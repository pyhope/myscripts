import argparse

parser = argparse.ArgumentParser(description="Extract a frame from LAMMPS dump file.")
parser.add_argument("--input_filename", "-i", type=str, default='npt.dump', help="Input filename (default: npt.dump)")
parser.add_argument("--output_filename", "-o", type=str, default='selected.dump', help="Output filename (default: selected.dump)")
parser.add_argument("--timestep", "-t", type=int, help="Timestep of the target frame. The last frame will be used if not provided.")

args = parser.parse_args()

with open(args.input_filename, 'r') as datain:
    last_frame = []
    target_frame = None
    line = datain.readline()

    while line:
        if "ITEM: TIMESTEP" in line:
            timestep_line = datain.readline()
            thisframe = int(timestep_line.strip())
            
            if target_frame is not None:
                break  # Exit if we have already found the target frame
            
            current_frame = ["ITEM: TIMESTEP\n", timestep_line]
            
            # Read the current frame
            while True:
                line = datain.readline()
                if not line or "ITEM: TIMESTEP" in line:
                    break
                current_frame.append(line)
                
            # Check if this is the target frame
            if args.timestep is not None and thisframe == args.timestep:
                target_frame = current_frame
            else:
                last_frame = current_frame
            
        else:
            line = datain.readline()

# Use the last frame if no specific timestep is provided or target frame is not found
if target_frame is None:
    target_frame = last_frame

if target_frame:
    with open(args.output_filename, 'w') as dataout:
        dataout.writelines(target_frame)
else:
    print(f"Timestep {args.timestep} not found in the input file.")
