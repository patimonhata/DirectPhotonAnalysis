#include "internal/TopologyEventRenderer.h"

#include <TFile.h>
#include <TSystem.h>

#include <iostream>
#include <string>

void BrowseTopologyEvents(
    const char* input_file,
    const char* current_event_pdf = "current_event.pdf",
    int topology_filter = -1,
    int pathway_filter = -1,
    int first_event = -1)
{
  TFile* file = TFile::Open(input_file, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "BrowseTopologyEvents - cannot open " << input_file << std::endl;
    return;
  }
  const auto events = topology_display::event_ids(file, topology_filter, pathway_filter);
  if (events.empty())
  {
    std::cerr << "BrowseTopologyEvents - no events match the filters" << std::endl;
    file->Close();
    return;
  }
  std::size_t position = 0;
  if (first_event >= 0)
  {
    const auto found = std::find(events.begin(), events.end(), first_event);
    if (found != events.end()) position = static_cast<std::size_t>(found - events.begin());
  }
  while (position < events.size())
  {
    const int event = events[position];
    if (!topology_display::print_event_pages(
            file, event, current_event_pdf, true, true))
    {
      std::cerr << "BrowseTopologyEvents - failed to render event " << event
                << std::endl;
      break;
    }
    gSystem->ProcessEvents();
    std::cout << "Rendered event " << event << " (" << position + 1 << "/"
              << events.size() << ") to " << current_event_pdf
              << ". Enter: next, p: previous, q: quit > " << std::flush;
    std::string command;
    std::getline(std::cin, command);
    if (command == "q" || command == "Q") break;
    if ((command == "p" || command == "P") && position > 0)
      --position;
    else
      ++position;
  }
  file->Close();
}
