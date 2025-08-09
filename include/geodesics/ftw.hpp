#include <cstdio>
#include <vector>
#include <charconv>
#include <fstream>

// one-shot buffered text writer (newline-terminated rows)
struct FastTextWriter {
  FILE* f = nullptr;
  std::vector<char> buf;
  size_t pos = 0;

  explicit FastTextWriter(const char* path, size_t buf_bytes = 8u << 20) : buf(buf_bytes) {
    f = std::fopen(path, "wb");
    if (f) setvbuf(f, nullptr, _IOFBF, 1 << 20);
  }
  ~FastTextWriter() { flush(); if (f) std::fclose(f); }

  void flush() {
    if (f && pos) { std::fwrite(buf.data(), 1, pos, f); pos = 0; }
  }
  inline void put(char c) { if (pos == buf.size()) flush(); buf[pos++] = c; }

  // print a single double with \n at the end
  inline void print_double(double x) {
    if (buf.size() - pos < 32) flush();
    auto r = std::to_chars(buf.data() + pos, buf.data() + pos + 32, x,
      std::chars_format::general, 17);
    pos = size_t(r.ptr - buf.data());
  }
  template<class... Ts>
  void row(Ts... xs) { ((print_double(xs), put(' ')), ...); buf[pos - 1] = '\n'; }
};