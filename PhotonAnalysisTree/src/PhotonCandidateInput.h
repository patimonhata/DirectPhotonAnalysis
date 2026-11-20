#ifndef PHOTON_TREE_PHOTON_CANDIDATE_INPUT_H
#define PHOTON_TREE_PHOTON_CANDIDATE_INPUT_H

#include <exception>
#include <limits>
#include <string>

namespace photon_tree
{
inline bool photon_candidate_source_segment(const std::string& filename, unsigned int& segment)
{
  const auto basename_begin = filename.find_last_of("/\\");
  const auto basename = filename.substr(basename_begin == std::string::npos ? 0 : basename_begin + 1);
  const auto extension = basename.rfind(".root");
  const auto separator = extension == std::string::npos ? std::string::npos : basename.rfind('-', extension);
  if (extension == std::string::npos || extension + 5 != basename.size() || separator == std::string::npos || separator + 1 >= extension)
    return false;

  const auto digits = basename.substr(separator + 1, extension - separator - 1);
  std::size_t parsed = 0;
  unsigned long value = 0;
  try
  {
    value = std::stoul(digits, &parsed);
  }
  catch (const std::exception&)
  {
    return false;
  }
  if (parsed != digits.size() || value > std::numeric_limits<unsigned int>::max()) return false;
  segment = static_cast<unsigned int>(value);
  return true;
}

inline unsigned long long photon_candidate_event_uid(unsigned int source_file_id, unsigned int event_in_file)
{
  return (static_cast<unsigned long long>(source_file_id) << 32U) | event_in_file;
}
}

#endif
