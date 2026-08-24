#include "internal/TopologyEventRenderer.h"

#include <TFile.h>

#include <algorithm>
#include <iostream>
#include <string>

void MakeTopologyEventBook(
    const char* input_file,
    const char* output_pdf = "topology_event_book.pdf",
    int topology_filter = -1,
    int pathway_filter = -1,
    int max_events = -1)
{
  TFile* file = TFile::Open(input_file, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "MakeTopologyEventBook - cannot open " << input_file << std::endl;
    return;
  }
  auto events = topology_display::event_ids(file, topology_filter, pathway_filter);
  if (max_events >= 0 && static_cast<int>(events.size()) > max_events)
    events.resize(static_cast<std::size_t>(max_events));
  if (events.empty())
  {
    std::cerr << "MakeTopologyEventBook - no events match the filters" << std::endl;
    file->Close();
    return;
  }
  for (std::size_t position = 0; position < events.size(); ++position)
  {
    const bool first = position == 0;
    const bool last = position + 1 == events.size();
    if (!topology_display::print_event_pages(
            file, events[position], output_pdf, first, last))
    {
      std::cerr << "MakeTopologyEventBook - failed at event "
                << events[position] << std::endl;
      file->Close();
      return;
    }
  }
  std::cout << "MakeTopologyEventBook - wrote " << events.size()
            << " events to " << output_pdf << std::endl;
  file->Close();
}
