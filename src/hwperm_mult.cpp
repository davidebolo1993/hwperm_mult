#include <algorithm>
#include <cctype>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Options {
  std::string input_path;
  std::size_t nperm = 20000;
  double eps = 1e-10;
  char sep = '\t';
  bool quiet = false;
  bool print_matrix = false;
  bool seed_given = false;
  std::uint64_t seed = 0;
};

[[noreturn]] void fail(const std::string& msg) {
  throw std::runtime_error(msg);
}

void usage(std::ostream& os) {
  os << "Usage: hwperm_mult --input <file> [options]\n"
        "       cat file.tsv | hwperm_mult_fast_cli [options]\n\n"
        "Options:\n"
        "  --input <path>       File with 2 lines (line1 and line2).\n"
        "  --nperm <int>        Number of permutations (default: 20000).\n"
        "  --eps <double>       Near-equality tolerance (default: 1e-10).\n"
        "  --sep <token>        Separator: tab|space|comma|semicolon|<char> (default: tab).\n"
        "  --seed <int>         Fixed RNG seed (optional).\n"
        "  --quiet              Print only final key-value lines.\n"
        "  --print-matrix       Print triangular genotype matrix.\n"
        "  --help               Show this help.\n";
}

std::string trim(const std::string& s) {
  std::size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])) != 0) {
    ++start;
  }
  std::size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
    --end;
  }
  return s.substr(start, end - start);
}

std::string lower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

char parse_sep(const std::string& token) {
  const std::string t = lower_copy(token);
  if (t == "tab" || token == "\\t") {
    return '\t';
  }
  if (t == "space") {
    return ' ';
  }
  if (t == "comma") {
    return ',';
  }
  if (t == "semicolon") {
    return ';';
  }
  if (token.size() == 1) {
    return token[0];
  }
  fail("Invalid --sep value: " + token);
}

Options parse_options(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--help" || a == "-h") {
      usage(std::cout);
      std::exit(0);
    }
    if (a == "--quiet") {
      opt.quiet = true;
      continue;
    }
    if (a == "--print-matrix") {
      opt.print_matrix = true;
      continue;
    }
    auto needs_value = [&](const std::string& key) {
      if (i + 1 >= argc) {
        fail("Missing value for " + key);
      }
      return std::string(argv[++i]);
    };

    if (a == "--input") {
      opt.input_path = needs_value(a);
    } else if (a == "--nperm") {
      const std::string v = needs_value(a);
      std::size_t pos = 0;
      const unsigned long long parsed = std::stoull(v, &pos);
      if (pos != v.size() || parsed == 0) {
        fail("Invalid --nperm value: " + v);
      }
      opt.nperm = static_cast<std::size_t>(parsed);
    } else if (a == "--eps") {
      const std::string v = needs_value(a);
      std::size_t pos = 0;
      opt.eps = std::stod(v, &pos);
      if (pos != v.size() || !(opt.eps > 0.0)) {
        fail("Invalid --eps value: " + v);
      }
    } else if (a == "--sep") {
      opt.sep = parse_sep(needs_value(a));
    } else if (a == "--seed") {
      const std::string v = needs_value(a);
      std::size_t pos = 0;
      opt.seed = static_cast<std::uint64_t>(std::stoull(v, &pos));
      if (pos != v.size()) {
        fail("Invalid --seed value: " + v);
      }
      opt.seed_given = true;
    } else {
      fail("Unknown argument: " + a + " (use --help)");
    }
  }
  return opt;
}

std::vector<std::string> split_line(const std::string& line, char sep) {
  std::vector<std::string> out;
  if (sep == ' ') {
    std::string cur;
    for (char c : line) {
      if (std::isspace(static_cast<unsigned char>(c)) != 0) {
        if (!cur.empty()) {
          out.push_back(cur);
          cur.clear();
        }
      } else {
        cur.push_back(c);
      }
    }
    if (!cur.empty()) {
      out.push_back(cur);
    }
    return out;
  }

  std::string cur;
  for (char c : line) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

bool is_missing_token(const std::string& token) {
  const std::string t = lower_copy(trim(token));
  return t.empty() || t == "na" || t == "nan" || t == ".";
}

inline int tri_index(int i, int j) {
  return (i * (i + 1)) / 2 + j;
}

bool nearly_equal_scalar(double a, double b, double eps) {
  if (a == b) {
    return true;
  }
  const double absA = std::fabs(a);
  const double absB = std::fabs(b);
  const double diff = std::fabs(a - b);
  if (a == 0.0 || b == 0.0 || (absA + absB < DBL_MIN)) {
    return diff < (eps * DBL_MIN);
  }
  const double den = std::min(absA + absB, static_cast<double>(DBL_MAX));
  return (diff / den) < eps;
}

struct TriangularData {
  std::vector<std::string> labels;
  std::vector<std::vector<int>> matrix;
};

TriangularData alleles_to_triangular(const std::vector<std::string>& line1,
                                     const std::vector<std::string>& line2) {
  if (line1.size() != line2.size()) {
    fail("line1 and line2 must have equal length.");
  }

  std::vector<std::pair<std::string, std::string>> pairs;
  pairs.reserve(line1.size());
  std::vector<std::string> all_labels;
  all_labels.reserve(line1.size() * 2);

  for (std::size_t i = 0; i < line1.size(); ++i) {
    if (is_missing_token(line1[i]) || is_missing_token(line2[i])) {
      continue;
    }
    const std::string a1 = "A" + trim(line1[i]);
    const std::string a2 = "A" + trim(line2[i]);
    pairs.emplace_back(a1, a2);
    all_labels.push_back(a1);
    all_labels.push_back(a2);
  }

  std::sort(all_labels.begin(), all_labels.end());
  all_labels.erase(std::unique(all_labels.begin(), all_labels.end()), all_labels.end());
  if (all_labels.empty()) {
    fail("No non-missing allele pairs found.");
  }

  std::unordered_map<std::string, int> idx;
  idx.reserve(all_labels.size());
  for (int i = 0; i < static_cast<int>(all_labels.size()); ++i) {
    idx[all_labels[static_cast<std::size_t>(i)]] = i;
  }

  const int k = static_cast<int>(all_labels.size());
  std::vector<std::vector<int>> M(static_cast<std::size_t>(k), std::vector<int>(static_cast<std::size_t>(k), 0));

  for (const auto& p : pairs) {
    int i = idx[p.first];
    int j = idx[p.second];
    if (i < j) {
      std::swap(i, j);
    }
    M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] += 1;
  }

  return TriangularData{all_labels, M};
}

struct PermResult {
  double pofthesample = std::numeric_limits<double>::quiet_NaN();
  double pval = std::numeric_limits<double>::quiet_NaN();
  std::size_t n_used = 0;
};

PermResult hw_perm_autosomal(const std::vector<std::vector<int>>& x,
                             std::size_t nperm,
                             double eps,
                             std::mt19937_64& rng) {
  const int k = static_cast<int>(x.size());
  if (k == 0 || static_cast<int>(x[0].size()) != k) {
    fail("Matrix must be square and non-empty.");
  }

  std::vector<int> row_sum(static_cast<std::size_t>(k), 0);
  std::vector<int> col_sum(static_cast<std::size_t>(k), 0);
  int n = 0;
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      const int v = x[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      if (v < 0) {
        fail("Genotype matrix cannot contain negative values.");
      }
      row_sum[static_cast<std::size_t>(i)] += v;
      col_sum[static_cast<std::size_t>(j)] += v;
      n += v;
    }
  }
  if (n <= 0) {
    fail("Genotype matrix has zero total count.");
  }

  std::vector<int> ac_counts(static_cast<std::size_t>(k), 0);
  int total_alleles = 0;
  for (int i = 0; i < k; ++i) {
    ac_counts[static_cast<std::size_t>(i)] = row_sum[static_cast<std::size_t>(i)] + col_sum[static_cast<std::size_t>(i)];
    total_alleles += ac_counts[static_cast<std::size_t>(i)];
  }
  if (total_alleles != 2 * n) {
    fail("Invalid genotype matrix: total allele count mismatch.");
  }

  double h_obs = 0.0;
  double gt_lgamma_obs = 0.0;
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j <= i; ++j) {
      const int v = x[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      gt_lgamma_obs += std::lgamma(static_cast<double>(v) + 1.0);
      if (i != j) {
        h_obs += static_cast<double>(v);
      }
    }
  }

  double constant_part = std::lgamma(static_cast<double>(n) + 1.0) -
                         std::lgamma(static_cast<double>(2 * n) + 1.0);
  for (int c : ac_counts) {
    constant_part += std::lgamma(static_cast<double>(c) + 1.0);
  }
  const double log2v = std::log(2.0);
  const double pofthesample = std::exp(constant_part + h_obs * log2v - gt_lgamma_obs);

  std::vector<double> lgamma_table(static_cast<std::size_t>(n + 1), 0.0);
  for (int i = 0; i <= n; ++i) {
    lgamma_table[static_cast<std::size_t>(i)] = std::lgamma(static_cast<double>(i) + 1.0);
  }

  std::vector<int> allele_pool;
  allele_pool.reserve(static_cast<std::size_t>(total_alleles));
  for (int i = 0; i < k; ++i) {
    allele_pool.insert(allele_pool.end(), ac_counts[static_cast<std::size_t>(i)], i);
  }
  std::vector<int> perm_pool = allele_pool;

  const int tri_len = (k * (k + 1)) / 2;
  std::vector<int> gt(static_cast<std::size_t>(tri_len), 0);
  std::vector<int> diag_idx(static_cast<std::size_t>(k), 0);
  for (int i = 0; i < k; ++i) {
    diag_idx[static_cast<std::size_t>(i)] = tri_index(i, i);
  }

  int nsmaller = 0;
  for (std::size_t rep = 0; rep < nperm; ++rep) {
    std::shuffle(perm_pool.begin(), perm_pool.end(), rng);
    std::fill(gt.begin(), gt.end(), 0);

    for (int i = 0; i < total_alleles; i += 2) {
      int a = perm_pool[static_cast<std::size_t>(i)];
      int b = perm_pool[static_cast<std::size_t>(i + 1)];
      if (a < b) {
        std::swap(a, b);
      }
      gt[static_cast<std::size_t>(tri_index(a, b))] += 1;
    }

    int homozyg = 0;
    double gt_lgamma = 0.0;
    for (int c : gt) {
      gt_lgamma += lgamma_table[static_cast<std::size_t>(c)];
    }
    for (int idx : diag_idx) {
      homozyg += gt[static_cast<std::size_t>(idx)];
    }

    const int heter = n - homozyg;
    const double d = std::exp(constant_part + static_cast<double>(heter) * log2v - gt_lgamma);
    if (nearly_equal_scalar(d, pofthesample, eps) || d < pofthesample) {
      ++nsmaller;
    }
  }

  const double pval = static_cast<double>(nsmaller) / static_cast<double>(nperm);
  return PermResult{pofthesample, pval, static_cast<std::size_t>(n)};
}

std::pair<std::string, std::string> read_two_lines(const std::string& path) {
  std::istream* in = &std::cin;
  std::ifstream f;
  if (!path.empty()) {
    f.open(path);
    if (!f.is_open()) {
      fail("Cannot open input file: " + path);
    }
    in = &f;
  }

  std::string line1;
  std::string line2;
  if (!std::getline(*in, line1) || !std::getline(*in, line2)) {
    fail("Expected at least 2 lines (line1 then line2).");
  }
  if (!line1.empty() && line1.back() == '\r') {
    line1.pop_back();
  }
  if (!line2.empty() && line2.back() == '\r') {
    line2.pop_back();
  }
  return std::make_pair(line1, line2);
}

} // namespace

int main(int argc, char** argv) {
  try {
    const Options opt = parse_options(argc, argv);
    const auto lines = read_two_lines(opt.input_path);
    const std::vector<std::string> line1 = split_line(lines.first, opt.sep);
    const std::vector<std::string> line2 = split_line(lines.second, opt.sep);

    TriangularData tri = alleles_to_triangular(line1, line2);

    std::mt19937_64 rng;
    if (opt.seed_given) {
      rng.seed(opt.seed);
    } else {
      std::random_device rd;
      std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
      rng.seed(seq);
    }

    const PermResult res = hw_perm_autosomal(tri.matrix, opt.nperm, opt.eps, rng);

    if (!opt.quiet) {
      std::cout << "Permutation test for Hardy-Weinberg equilibrium (autosomal)." << '\n';
      std::cout << tri.labels.size() << " alleles detected." << '\n';
      std::cout << "Observed statistic: " << std::setprecision(12) << res.pofthesample
                << "  " << opt.nperm << " permutations. p-value: " << res.pval << '\n';
    }

    std::cout << "p-value: " << std::setprecision(12) << res.pval << '\n';
    std::cout << "observed_statistic: " << std::setprecision(12) << res.pofthesample << '\n';
    std::cout << "k_alleles: " << tri.labels.size() << '\n';
    std::cout << "n_genotypes: " << res.n_used << '\n';

    if (opt.print_matrix) {
      std::cout << "genotype_matrix_lower_triangular:" << '\n';
      for (std::size_t i = 0; i < tri.labels.size(); ++i) {
        for (std::size_t j = 0; j < tri.labels.size(); ++j) {
          if (j != 0) {
            std::cout << '\t';
          }
          std::cout << tri.matrix[i][j];
        }
        std::cout << '\n';
      }
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
}
