#ifndef ETA_ORIGIN_FRACTIONS_H
#define ETA_ORIGIN_FRACTIONS_H
#include <cmath>
#include <map>
#include <utility>
#include <vector>
namespace eta_origin
{
struct Majority
{
  double fraction = 0.0, pi0_total_fraction = 0.0;
  int mode = 0, barcode = 0, track = 0;
  bool unresolved = false;
};
inline Majority summarize(unsigned int begin, unsigned int end, const std::vector<int>& embedding, const std::vector<float>& fraction,
                          const std::vector<int>& barcode, const std::vector<int>& track, const std::vector<int>& mode,
                          const std::vector<int>& split_status, const std::vector<int>& pi0_barcode, const std::vector<int>& pi0_track)
{
  struct Group { double fraction = 0; int mode = 0; bool unresolved = false; };
  std::map<std::pair<int, int>, Group> groups;
  Majority result;
  for (unsigned int i = begin; i < end; ++i)
  {
    if (embedding[i] != 1) continue;
    if (pi0_barcode[i] != 0 || pi0_track[i] != 0) result.pi0_total_fraction += fraction[i];
    if (barcode[i] == 0 && track[i] == 0) continue;
    // A G4-primary eta and its HepMC identity belong to one parent, not two groups.
    const auto key = barcode[i] != 0 ? std::make_pair(barcode[i], 0) : std::make_pair(0, track[i]);
    auto& group = groups[key];
    if (group.fraction > 0 && group.mode != mode[i]) group.unresolved = true;
    group.mode = mode[i]; group.unresolved = group.unresolved || split_status[i] < 0;
    group.fraction += fraction[i];
  }
  for (const auto& [key, group] : groups) if (group.fraction > result.fraction)
  { result.fraction = group.fraction; result.mode = group.mode; result.barcode = key.first; result.track = key.second; result.unresolved = group.unresolved; }
  return result;
}
}
#endif
