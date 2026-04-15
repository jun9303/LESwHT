#! /usr/bin/env python3
import argparse
import os
from dimple_generator import (
    DimplePatternGenerator,
    case_name_from_params,
    prepare_case_directories,
    write_case_env,
)

def main():
    # Set up the argument parser with RawTextHelpFormatter to allow multi-line descriptions
    parser = argparse.ArgumentParser(
        description="Generate repetitive dimple surface patterns from baseline STL files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Define CLI arguments
    parser.add_argument(
        '-t', '--topology', 
        type=int, 
        choices=range(7), 
        default=0,
        help=(
            "Baseline topology ID [0-6]:\n"
            "  0 : Cylinder (Circular top, cylindrical carve)\n"
            "  1 : Spherical (Circular top, smooth depth variation)\n"
            "  2 : Teardrop Down (Isosceles + circle shape)\n"
            "  3 : Teardrop Up (X-flipped Teardrop Down)\n"
            "  4 : Diamond (Overlapping teardrops)\n"
            "  5 : Triangle Down (Isosceles shape, linear depth)\n"
            "  6 : Triangle Up (X-flipped Triangle Down)"
        )
    )
    
    parser.add_argument(
        '-d', '--depth', 
        type=float, 
        default=30.0,
        help="Dimple depth in wall units (e.g., 10, 15, 20, 25, 30, 35, 40, 45, 50)"
    )
    
    parser.add_argument(
        '-s', '--scale', 
        type=float, 
        default=1.0,
        help="Isotropic pattern scaling factor (e.g., 0.2, 0.4, 0.6, 0.8, 1.0)\n(1.0 means neighboring dimple tops are exactly tangential)"
    )
    
    parser.add_argument(
        '-x', '--stream_scale', 
        type=float, 
        default=1.0,
        help=(
            "Streamwise scaling factor (e.g., 0.5, 0.75, 1.0, 1.25, 1.5). "
            "Stretches each dimple in x and scales streamwise pattern spacing."
        )
    )

    parser.add_argument(
        '--compact',
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Generate a compact periodic cell: z in [-5, 5] and x spanning 2 upstream/2 downstream "
            "pattern spacings around the center."
        )
    )

    parser.add_argument(
        '--output-root',
        type=str,
        default='../output',
        help='Root folder where per-run case directories are created.'
    )

    args = parser.parse_args()

    # Initialize and run the generator
    generator = DimplePatternGenerator("baselineSTL")
    
    try:
        case_name = case_name_from_params(args.topology, args.depth, args.scale, args.stream_scale)
        case_root = os.path.abspath(os.path.join(args.output_root, case_name))
        prepare_case_directories(case_root)

        stl_path = os.path.join(case_root, 'geometry', 'dimple_slab_ascii.stl')
        generator.generate_to_path(
            args.topology,
            args.depth,
            args.scale,
            args.stream_scale,
            stl_path,
            compact=args.compact,
        )

        write_case_env(case_root, case_name, stl_path)
        print(f"Case output directory prepared: {case_root}")
        print(f"Run context written to: {os.path.abspath(os.path.join(args.output_root, 'running_case.env'))}")
    except FileNotFoundError as e:
        print(f"\nError: Could not find baseline STL file. {e}")
        print("Make sure your baseline files are correctly named in the 'baselineSTL' directory.\n")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}\n")

if __name__ == "__main__":
    main()