#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

parser = argparse.ArgumentParser(description="Calculate the average distance between helium atoms")
parser.add_argument("--input_file", "-i", type=str, default="./XDATCAR",  help="Path to the XDATCAR file")
parser.add_argument("--frequency", "-f", type=int, default=100, help="Frequency of configurations to extract")
args = parser.parse_args()

def extract_configurations(file_content):
    # Split the content by the configuration header
    splits = file_content.split("Direct configuration=")
    head = splits[0]
    configurations = splits[1:]

    # Filter configurations every 200
    desired_configurations = {i: config for i, config in enumerate(configurations, 1) if i % args.frequency == 0}
    
    return head, desired_configurations

def write_to_files(head, configurations):
    parent_dir = "Extracted_Configurations"
    
    # Ensure the directory exists or create it
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)
    
    for i, config in configurations.items():
        directory_name = f"{i}"
        directory_path = os.path.join(parent_dir, directory_name)
        
        # Ensure the directory exists or create it
        if not os.path.exists(directory_path):
            os.mkdir(directory_path)
        
        with open(os.path.join(directory_path, "POSCAR"), 'w') as f:
            # Write the head and the current configuration
            f.write(head)
            f.write("Direct configuration=" + config)

if __name__ == "__main__":
    with open(args.input_file, 'r') as f:
        content = f.read()

    head, configs = extract_configurations(content)
    write_to_files(head, configs)
