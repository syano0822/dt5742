import ROOT
import os
import sys

def merge_and_compare(file0_path, file1_path, output_dir, num_events=5):
    """
    Reads two waveform_analyzed.root files (with 'Analysis' TTree),
    extracts ampMax for matching events, and produces merged Heatmap plots.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Opening File 0: {file0_path}")
    f0 = ROOT.TFile.Open(file0_path, "READ")
    print(f"Opening File 1: {file1_path}")
    f1 = ROOT.TFile.Open(file1_path, "READ")

    if not f0 or f0.IsZombie():
        print(f"Error opening {file0_path}")
        return
    if not f1 or f1.IsZombie():
        print(f"Error opening {file1_path}")
        return

    # Get Trees
    # Default tree name from config is "Analysis"
    tree0 = f0.Get("Analysis")
    tree1 = f1.Get("Analysis")

    if not tree0:
        print("Error: 'Analysis' tree not found in file 0")
        return
    if not tree1:
        print("Error: 'Analysis' tree not found in file 1")
        return

    entries0 = tree0.GetEntries()
    entries1 = tree1.GetEntries()
    print(f"Entries: File0={entries0}, File1={entries1}")

    # Process first N events
    # Assuming events are synchronized and sorted (1-to-1 mapping by index)
    # If not, we would need to build an index. For QA, we assume 1-to-1 for now.
    
    max_process = min(entries0, entries1, num_events)
    print(f"Processing first {max_process} events...")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)

    for i in range(max_process):
        tree0.GetEntry(i)
        tree1.GetEntry(i)

        evt0 = getattr(tree0, "event", -1)
        evt1 = getattr(tree1, "event", -1)

        if evt0 != evt1:
            print(f"Warning: Event mismatch at index {i}: {evt0} vs {evt1}. Skipping.")
            continue
        
        print(f"Processing Event {evt0}...")

        # Collect data points for this event from both files
        # Structure: sensor_id -> list of (col, strip, amp)
        sensor_data = {}

        def extract_data(tree):
            # Vectors: sensorID, sensorCol (Strip), sensorRow (Col), ampMax, isHorizontal
            try:
                ids = tree.sensorID
                strips = tree.sensorCol
                cols = tree.sensorRow
                amps = tree.ampMax
                horiz = tree.isHorizontal # Use checking 0th element if needed
                
                n_ch = len(ids)
                for j in range(n_ch):
                    sid = ids[j]
                    s = strips[j]
                    c = cols[j]
                    a = amps[j]
                    
                    if sid not in sensor_data:
                        sensor_data[sid] = {'points': [], 'is_horiz': False} # default vert
                    
                    sensor_data[sid]['points'].append( (c, s, a) )
                    if len(horiz) > j:
                        sensor_data[sid]['is_horiz'] = bool(horiz[j])

            except Exception as e:
                print(f"Error extracting data from tree: {e}")

        extract_data(tree0)
        extract_data(tree1)

        # Now generate plots for each sensor
        sorted_sids = sorted(sensor_data.keys())
        num_sensors = len(sorted_sids)
        
        if num_sensors == 0:
            print(f"  No sensor data found for event {evt0}.")
            continue

        c_out = ROOT.TCanvas(f"merged_event_{evt0}", f"Merged Event {evt0}", 600 * num_sensors, 600)
        c_out.Divide(num_sensors, 1)

        histograms = [] # Keep objects alive

        for pad_idx, sid in enumerate(sorted_sids):
            c_out.cd(pad_idx + 1)
            
            data = sensor_data[sid]
            points = data['points']
            is_horiz = data['is_horiz']

            if not points:
                continue

            # Determine bounds
            # Vertical: X=Column, Y=Strip
            # Horizontal: X=Strip, Y=Column
            
            # Extract raw coords to find min/max
            raw_cols = [p[0] for p in points]
            raw_strips = [p[1] for p in points]
            
            if not is_horiz:
                # X = Col, Y = Strip
                x_vals = raw_cols
                y_vals = raw_strips
                x_title = "Column"
                y_title = "Strip"
            else:
                # X = Strip, Y = Col
                x_vals = raw_strips
                y_vals = raw_cols
                x_title = "Strip"
                y_title = "Column"

            min_x, max_x = min(x_vals), max(x_vals)
            min_y, max_y = min(y_vals), max(y_vals)

            # Create Hist
            h_name = f"hist_sensor{sid}_ev{evt0}"
            h = ROOT.TH2F(h_name, 
                          f"Sensor {sid} - Event {evt0};{x_title};{y_title};Amplitude (ADC)",
                          max_x - min_x + 1, min_x - 0.5, max_x + 0.5,
                          max_y - min_y + 1, min_y - 0.5, max_y + 0.5)
            
            # Fill
            for p in points:
                # p = (col, strip, amp)
                col, strip, amp = p
                if not is_horiz:
                    h.Fill(col, strip, amp)
                else:
                    h.Fill(strip, col, amp)
            
            h.SetMinimum(0)
            h.SetMaximum(300) # Fixed scale for consistency
            h.Draw("COLZ TEXT")
            histograms.append(h)

        out_name = os.path.join(output_dir, f"merged_event_{evt0:06d}_quality_check.png")
        c_out.SaveAs(out_name)
        print(f"  Saved {out_name}")

    f0.Close()
    f1.Close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 merge_qa_plots.py <daq00_analyzed.root> <daq01_analyzed.root> <output_dir>")
        sys.exit(1)
    
    merge_and_compare(sys.argv[1], sys.argv[2], sys.argv[3])