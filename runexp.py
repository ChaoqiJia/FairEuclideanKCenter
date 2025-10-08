import os
from time import time as time
import numpy as np
import csv
import math
from scipy.spatial import distance
from ouralg import one_two_sqrt_three_m_Approx
from collections import Counter
from copy import deepcopy
import re
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

numIter = 10  # Number of iterations
percentEachGroup = 20  # Percentage of members for each group
base = 0.0001
dataset = ["twitch"]
approx = 0.5


def norm():
    df = pd.read_csv(PATH_o, header=None)

    features = df.iloc[:, :-1] 
    labels = df.iloc[:, -1]  
    scaler = MinMaxScaler()
    normalized_features = scaler.fit_transform(features)
    normalized_df = pd.DataFrame(normalized_features, columns=features.columns)
    normalized_df[''] = labels  
    print(normalized_df.head())  
    normalized_df.to_csv(PATH, index=False, header=False)  


def evaluate(X, C):
    min_distances = np.min(distance.cdist(X, np.array(C), metric="euclidean"), axis=1)
    # min_distances = np.min(distance.cdist(X, ), axis=1)
    return max(min_distances)

def extract_constraints(file_path, target_filename):
    constraints_list = []
    pattern = re.compile(
        r'path:\s*(?P<path>[^;]+);\s*constraints:\s*\[(?P<constraints>[^\]]+)\]'
    )

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                path = match.group('path').strip()
                constraints_str = match.group('constraints').strip()
                if target_filename.__contains__(path):
                    constraints = [int(num) for num in constraints_str.split()]
                    break
                    # constraints_list.extend(constraints)

    return constraints


def extract_constraints_real(file_path):
    df = pd.read_csv(file_path)  
    last_column = df.iloc[:, -1]
    counts = last_column.value_counts()
    return counts


############## real-world dataset/ approximation ratio #####################
def Simu_approx(iterN):
    epsilon = 0.05
    for r in range(len(dataset)):
        for group_n in range(1):
            for i_iter in range(iterN):
                PATH = "Dataset/" + dataset[r] + f"_norm.csv"
                PATH_G = "Dataset/" + dataset[r] + f"_group.csv"
                # counts_m = extract_constraints(PATH_oo_constrain_file, PATH)
                counts_m = extract_constraints_real(PATH_G)
                totalTime = []
                totalLoss = []
                constraints = [math.ceil(x / sum(counts_m) * percentEachGroup) for x in counts_m] 
                # constraints = counts_m
                print(constraints)

                start = time()
                C =  one_two_sqrt_three_m_Approx(constraints, epsilon, PATH, PATH_G)
                
                totalTime.append(time() - start)
                max_cost = 0
                
                with open(PATH, 'r') as csvfile:
                    csv_reader = csv.reader(csvfile)
                    next(csv_reader)
                    for row in csv_reader:
                        X = [float(item) for item in row[:-1]]
                        cost = evaluate([X], C)
                        max_cost = max(max_cost, cost)
                
                totalLoss.append(max_cost / approx)


                data = {
                    'dataset': dataset[r],
                    'alg': ["1+2√3-Approx"],
                    'approx': totalLoss,
                    'runtime': totalTime,
                    'm': len(constraints),
                    'k': [constraints]
                }
                print(data)

                df = pd.DataFrame(data)

                file_exists = os.path.isfile('output/test.csv')
                df.to_csv('output/test.csv', mode='a', index=False, header=not file_exists)

Simu_approx(numIter)
