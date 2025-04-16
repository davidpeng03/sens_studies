import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class LHEParser:
    def __init__(self, lhe_file):
        self.lhe_file = lhe_file
        self.events_data = []

    def parse_lhe(self):
        """
        Parse the LHE file and extract events into a DataFrame.
        """
        print(f"Parsing LHE file: {self.lhe_file}...")
        
        with open(self.lhe_file, "r") as file:
            inside_event = False
            event_lines = []

            for line in file:
                line = line.strip()
                if line.startswith("<event>"):
                    inside_event = True
                    event_lines = []
                elif line.startswith("</event>"):
                    inside_event = False
                    self._process_event(event_lines)
                elif inside_event:
                    event_lines.append(line)
        
        print(f"Extracted {len(self.events_data)} events.")
        return self.to_dataframe()

    def _process_event(self, event_lines):
        """
        Process a single event and extract relevant data.
        """
        if not event_lines:
            return

        header_parts = event_lines[0].split()
        num_particles = int(header_parts[0])
        process_id = int(header_parts[1])
        event_weight = float(header_parts[2])
        scale_Q = float(header_parts[3])
        alpha_em = float(header_parts[4])
        alpha_s = float(header_parts[5])

        event_info = []
        mothers_map = {}
        daughters_map = {}
        
        for i in range(1, num_particles + 1):
            parts = event_lines[i].split()
            pid = int(parts[0])
            status = int(parts[1])
            mother1 = int(parts[2])
            mother2 = int(parts[3])
            color1 = int(parts[4])
            color2 = int(parts[5])
            px, py, pz, E = map(float, parts[6:10])
            mass = float(parts[10])
            lifetime = float(parts[11])
            helicity = float(parts[12])
            pt = np.sqrt(px**2 + py**2)
            eta = 0.5 * np.log((E + pz) / (E - pz)) if (E - pz) != 0 else np.nan
            phi = np.arctan2(py, px)

            if mother1 > 0:
                mothers_map.setdefault(mother1, []).append(pid)
            if mother2 > 0:
                mothers_map.setdefault(mother2, []).append(pid)

            daughters_map[pid] = (mother1, mother2)

            event_info.append({
                "event_id": len(self.events_data) + 1,
                "particle_id": pid,
                "status": status,
                "mother1": mother1,
                "mother2": mother2,
                "color1": color1,
                "color2": color2,
                "px": px,
                "py": py,
                "pz": pz,
                "E": E,
                "mass": mass,
                "pT": pt,
                "eta": eta,
                "phi": phi,
                "lifetime": lifetime,
                "helicity": helicity,
                "process_id": process_id,
                "event_weight": event_weight,
                "scale_Q": scale_Q,
                "alpha_em": alpha_em,
                "alpha_s": alpha_s,
                "daughters": mothers_map.get(pid, []),
                "mothers": daughters_map.get(pid, (None, None))
            })
        
        self.events_data.append(event_info)

    def to_dataframe(self):
        """
        Convert parsed events to a Pandas DataFrame.
        """
        flat_data = [item for event in self.events_data for item in event]
        df = pd.DataFrame(flat_data)
        return df


if __name__ == "__main__":
    lhe_file = "db/Temp/madgraph/Events/run_02/unweighted_events.lhe" 
    parser = LHEParser(lhe_file)
    df = parser.parse_lhe()
    print(df.head(10))
    df.to_csv("df.csv")
    particle_id_of_interest = 9900012
    df_filtered = df[df["particle_id"] == particle_id_of_interest]
    
    plot_features = ["E", "pT", "pz", "eta", "phi", "mass", "event_weight", "scale_Q", "alpha_em", "alpha_s"]
    
    for feature in plot_features:
        plt.figure()
        plt.hist(df_filtered[feature], bins=30, alpha=0.75)
        plt.xlabel(feature)
        plt.ylabel("Count")
        plt.title(f"Distribution of {feature} for particle {particle_id_of_interest}")
        plt.show()
    
    plt.figure()
    df_filtered["mother1"].value_counts().plot(kind='bar', alpha=0.75)
    plt.xlabel("Mother Particle ID")
    plt.ylabel("Count")
    plt.title(f"Mother particles of {particle_id_of_interest}")
    plt.show()
    
    # daughters_list = [daughter for sublist in df_filtered["daughters"] for daughter in sublist]
    # plt.figure()
    # pd.Series(daughters_list).value_counts().plot(kind='bar', alpha=0.75)
    # plt.xlabel("Daughter Particle ID")
    # plt.ylabel("Count")
    # plt.title(f"Répartition des particules filles du {particle_id_of_interest}")
    # plt.show()

