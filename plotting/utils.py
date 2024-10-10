import numpy as np

def read_probs_data(file_path):
    """
    This function reads RNA BPPs or alignment coinc probs from a file.
    
    Args:
    file_path (str): Path to the file containing base pair probabilities.
    
    Returns:
    data (list of tuples): Each tuple contains (i, j, prob)
    """
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                i, j, prob = int(parts[0]), int(parts[1]), float(parts[2])
                data.append((i, j, prob if prob > 0.1 else np.exp(2 * prob) - 1))
    return data
