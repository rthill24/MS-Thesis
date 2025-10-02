import subprocess
import os

# List of (script_path, script_dir) tuples
scripts = [
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 1/SD_Midship_Section_Optimizer_Case_1.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 1"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 2/SD_Midship_Section_Optimizer_Case_2.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 2"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 3/SD_Midship_Section_Optimizer_Case_3.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 3"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 4/SD_Midship_Section_Optimizer_Case_4.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 4"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 1/SD_Midship_Section_Post_Processing_Case_1.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 1"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 2/SD_Midship_Section_Post_Processing_Case_2.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 2"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 3/SD_Midship_Section_Post_Processing_Case_3.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 3"),
    ("C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 4/SD_Midship_Section_Post_Processing_Case_4.py", "C:/Users/rthill/Documents/MS-Thesis/Midship Section Optimization/RB & Return to Optimality/Case 4"),
]

for script_path, script_dir in scripts:
    # Output file path in script's directory
    out_file = os.path.join(script_dir, "output.txt")
    with open(out_file, "w") as out:
        # Set cwd so the script runs in its own directory
        subprocess.run(
        ["python", os.path.basename(script_path)],
        cwd=script_dir
)