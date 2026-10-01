#include "NESTObservableArray.hh"

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace {
// Call f once per field, passing that field of each array. Every per-event
// field must be listed here; the methods below all go through this list.
template <typename F, typename... Arrays>
void for_each_field(F&& f, Arrays&&... arrays) {
  f(std::forward<Arrays>(arrays).energy_kev...);
  f(std::forward<Arrays>(arrays).x_mm...);
  f(std::forward<Arrays>(arrays).y_mm...);
  f(std::forward<Arrays>(arrays).z_mm...);
  f(std::forward<Arrays>(arrays).n_electrons...);
  f(std::forward<Arrays>(arrays).n_photons...);
  f(std::forward<Arrays>(arrays).s1_nhits...);
  f(std::forward<Arrays>(arrays).s1_nhits_dpe...);
  f(std::forward<Arrays>(arrays).s1_nhits_thr...);
  f(std::forward<Arrays>(arrays).s1r_phe...);
  f(std::forward<Arrays>(arrays).s1c_phe...);
  f(std::forward<Arrays>(arrays).s1r_phd...);
  f(std::forward<Arrays>(arrays).s1c_phd...);
  f(std::forward<Arrays>(arrays).s1r_spike...);
  f(std::forward<Arrays>(arrays).s1c_spike...);
  f(std::forward<Arrays>(arrays).s1_photon_times...);
  f(std::forward<Arrays>(arrays).s2_Nee...);
  f(std::forward<Arrays>(arrays).s2_Nph...);
  f(std::forward<Arrays>(arrays).s2_nhits...);
  f(std::forward<Arrays>(arrays).s2_nhits_dpe...);
  f(std::forward<Arrays>(arrays).s2r_phe...);
  f(std::forward<Arrays>(arrays).s2c_phe...);
  f(std::forward<Arrays>(arrays).s2r_phd...);
  f(std::forward<Arrays>(arrays).s2c_phd...);
  f(std::forward<Arrays>(arrays).s1_waveform_time...);
  f(std::forward<Arrays>(arrays).s1_waveform_amp...);
  f(std::forward<Arrays>(arrays).s2_waveform_time...);
  f(std::forward<Arrays>(arrays).s2_waveform_amp...);
}

template <typename T>
void append_field(std::vector<T>& dst, const std::vector<T>& src) {
  dst.insert(dst.end(), src.begin(), src.end());
}

template <typename T>
void append_field(std::vector<T>& dst, std::vector<T>&& src) {
  if (dst.empty())
    dst = std::move(src);
  else
    dst.insert(dst.end(), std::make_move_iterator(src.begin()),
               std::make_move_iterator(src.end()));
}

// Src is a const& (fields are copied) or an rvalue (fields are moved)
template <typename Src>
void append_all_fields(NESTObservableArray& dst, Src&& src) {
  for_each_field(
      [](auto& d, auto&& s) { append_field(d, std::forward<decltype(s)>(s)); },
      dst, std::forward<Src>(src));
}
}  // namespace

NESTObservableArray& NESTObservableArray::operator+=(const NESTObservableArray& other) {
  if (&other == this) {
    // Inserting a vector's own range into itself is undefined, so copy first
    NESTObservableArray copy(other);
    return *this += std::move(copy);
  }
  append_all_fields(*this, other);
  return *this;
}

NESTObservableArray& NESTObservableArray::operator+=(NESTObservableArray&& other) {
  if (&other == this) return *this += static_cast<const NESTObservableArray&>(other);
  append_all_fields(*this, std::move(other));
  return *this;
}

void NESTObservableArray::resize(std::size_t n_events) {
  for_each_field([n_events](auto& field) { field.resize(n_events); }, *this);
}

void NESTObservableArray::assign_events(std::size_t offset, NESTObservableArray&& src) {
  if (offset > size() or src.size() > size() - offset)
    throw std::out_of_range("NESTObservableArray::assign_events: events do not fit");
  // Only elements are written, never the vectors themselves, so concurrent
  // calls on disjoint ranges are safe
  for_each_field(
      [offset](auto& d, auto& s) { std::move(s.begin(), s.end(), d.begin() + offset); },
      *this, src);
}

// Method to store the results from simulating an S1 and an S2
void NESTObservableArray::store_signals(
  double energy,
  std::vector<double> pos,
  NEST::NESTresult result,
  std::vector<double> s1,
  std::vector<double> s2,
  std::vector<int64_t> s1_wf_time,
  std::vector<double> s1_wf_amp,
  std::vector<int64_t> s2_wf_time,
  std::vector<double> s2_wf_amp
){
  // Handle truth
  energy_kev.push_back(energy);
  x_mm.push_back(pos[0]);
  y_mm.push_back(pos[1]);
  z_mm.push_back(pos[2]);

  // Handle quanta
  n_photons.push_back(result.quanta.photons);
  n_electrons.push_back(result.quanta.electrons);
  s1_photon_times.push_back(result.photon_times);

  // Handle S1
  // 
  // MC-true integer hits in same OR different PMTs, NO double phe effect
  s1_nhits.push_back(std::abs(int(s1[0]))); 
  // MC-true integer hits WITH double phe effect (Nphe > nHits)
  s1_nhits_dpe.push_back(std::abs(int(s1[1])));
  // smeared DAQ pulse areas in phe, NO XYZ correction
  s1r_phe.push_back(std::abs(s1[2])); 
  // same as s1r_phe WITH XYZ correction
  s1c_phe.push_back(std::abs(s1[3])); 
  //  same as s1r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
  s1r_phd.push_back(std::abs(s1[4])); 
  // same as s1r_phd, but XYZ-corrected
  s1c_phd.push_back(std::abs(s1[5])); 
  // spike count, NO XYZ correction
  s1r_spike.push_back(std::abs(s1[6])); 
  // spike count, WITH XYZ correction
  s1c_spike.push_back(std::abs(s1[7]));   
  // nHits post coincidence window and single PE eff!
  s1_nhits_thr.push_back(std::abs(int(s1[8]))); 

  // Handle S2
  // 
  // integer number of electrons unabsorbed in liquid then getting extracted
  s2_Nee.push_back(std::abs(int(s2[0])));
  // raw number of photons produced in the gas gap
  s2_Nph.push_back(std::abs(int(s2[1])));  
  // MC-true integer hits in same OR different PMTs, NO double phe effect
  s2_nhits.push_back(std::abs(int(s2[2]))); 
  // MC-true integer hits WITH double phe effect (Nphe > nHits)
  s2_nhits_dpe.push_back(std::abs(int(s2[3]))); 
  // smeared DAQ pulse areas in phe, NO XYZ correction
  s2r_phe.push_back(std::abs(s2[4])); 
  // same as s2r_phe WITH XYZ correction
  s2c_phe.push_back(std::abs(s2[5])); 
  // same as s2r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
  s2r_phd.push_back(std::abs(s2[6])); 
  // same as s2r_phd, but XYZ-corrected
  s2c_phd.push_back(std::abs(s2[7])); 

  // Handle waveforms
  s1_waveform_time.push_back(s1_wf_time);
  s1_waveform_amp.push_back(s1_wf_amp);
  s2_waveform_time.push_back(s2_wf_time);
  s2_waveform_amp.push_back(s2_wf_amp);
}
