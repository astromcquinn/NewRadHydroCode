// A snapshot write that fails after its file was opened (disk full, quota
// exceeded) is reported and counted, and the run goes on
// (report_output_write_failure, io.hpp). /dev/full opens but fails every
// write with ENOSPC, the quota-full case of 2026-11-11 (empty snapshot and
// movie files written without a word). A good write is not counted.
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "check.hpp"
#include "radray/io.hpp"

using namespace radray;

int main() {
  const std::vector<float> v(1 << 20, 1.0f);  // 4 MB: more than any stream buffer
  const std::string tmp =
      (std::filesystem::temp_directory_path() / "radray_test_output_write.bin").string();
  write_raw_f32(tmp, v);
  CHECK(output_write_failures() == 0);
  CHECK(std::filesystem::file_size(tmp) == v.size() * sizeof(float));
  std::filesystem::remove(tmp);
  if (!std::filesystem::exists("/dev/full")) {
    std::printf("no /dev/full: the failure case is skipped\n");
    TEST_MAIN_END();
  }
  bool threw = false;
  try {
    write_raw_f32("/dev/full", v);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(!threw);
  CHECK(output_write_failures() == 1);
  write_raw_f32("/dev/full", v);
  CHECK(output_write_failures() == 2);
  TEST_MAIN_END();
}
