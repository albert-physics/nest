#ifndef NESTOBSERVABLEARRAY_HH
#define NESTOBSERVABLEARRAY_HH 1

#include <cstdint>
#include <utility>
#include <vector>

#include "NEST.hh"

// Per-event NEST outputs (truth, quanta, S1, S2, waveforms) stored as one vector per quantity
class NESTObservableArray {
  public:
    NESTObservableArray(){};
    ~NESTObservableArray() = default;

    void store_signals(
      double energy,
      std::vector<double> pos,
      NEST::NESTresult result,
      std::vector<double> s1,
      std::vector<double> s2,
      std::vector<int64_t> s1_wf_time,
      std::vector<double> s1_wf_amp,
      std::vector<int64_t> s2_wf_time,
      std::vector<double> s2_wf_amp
    );

    // Append the events of another array after this one's
    NESTObservableArray& operator+=(const NESTObservableArray& other);
    NESTObservableArray& operator+=(NESTObservableArray&& other);

    // Number of events stored
    std::size_t size() const { return energy_kev.size(); }

    // Resize every field to hold n_events (new events are value-initialised)
    void resize(std::size_t n_events);

    // Move the events of src into events [offset, offset + src.size()), which
    // must already exist (see resize). Calls writing disjoint ranges may run
    // concurrently from different threads.
    void assign_events(std::size_t offset, NESTObservableArray&& src);

    // Truth
    std::vector<double> energy_kev;
    std::vector<double> x_mm;
    std::vector<double> y_mm;
    std::vector<double> z_mm;

    // Quanta
    std::vector<int> n_electrons;
    std::vector<int> n_photons;

    // S1
    std::vector<int> s1_nhits; // MC-true integer hits in same OR different PMTs, NO double phe effect
    std::vector<int> s1_nhits_dpe; // MC-true integer hits WITH double phe effect (Nphe > nHits) 
    std::vector<int> s1_nhits_thr; // nHits post coincidence window and single PE eff! 
    std::vector<double> s1r_phe; // smeared DAQ pulse areas in phe, NO XYZ correction   
    std::vector<double> s1c_phe; // same as s1r_phe WITH XYZ correction 
    std::vector<double> s1r_phd; // same as s1r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)   
    std::vector<double> s1c_phd; // same as s1r_phd, but XYZ-corrected   
    std::vector<double> s1r_spike; // spike count, NO XYZ correction   
    std::vector<double> s1c_spike; // spike count, WITH XYZ correction   

    std::vector<std::vector<double>> s1_photon_times; // Time that S1 photons reach a PMT (includes undetected photons)

    // S2
    std::vector<int> s2_Nee; // integer number of electrons unabsorbed in liquid then getting extracted
    std::vector<int> s2_Nph; // raw number of photons produced in the gas gap
    std::vector<int> s2_nhits; // MC-true integer hits in same OR different PMTs, NO double phe effect
    std::vector<int> s2_nhits_dpe; // MC-true integer hits WITH double phe effect (Nphe > nHits)
    std::vector<double> s2r_phe; // smeared DAQ pulse areas in phe, NO XYZ correction
    std::vector<double> s2c_phe; // same as s2r_phe WITH XYZ correction
    std::vector<double> s2r_phd; // same s2r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
    std::vector<double> s2c_phd; // same as s2r_phd, but XYZ-corrected

    // Waveforms
    std::vector<std::vector<int64_t>> s1_waveform_time;
    std::vector<std::vector<double>> s1_waveform_amp;
    std::vector<std::vector<int64_t>> s2_waveform_time;
    std::vector<std::vector<double>> s2_waveform_amp;
};

// Concatenate two arrays: the events of rhs follow those of lhs
inline NESTObservableArray operator+(NESTObservableArray lhs, const NESTObservableArray& rhs) {
  lhs += rhs;
  return lhs;
}

inline NESTObservableArray operator+(NESTObservableArray lhs, NESTObservableArray&& rhs) {
  lhs += std::move(rhs);
  return lhs;
}

#endif
