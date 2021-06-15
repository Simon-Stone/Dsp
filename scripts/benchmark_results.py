import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

if __name__ == "__main__":
    res = pd.read_csv("../test/dsp_benchmark.csv", delimiter="\t")
    sns.boxplot(x=res['Task'], y=res['Execution time [us]'], hue=res['Implementation'])
    
    
    
    plt.show()
    