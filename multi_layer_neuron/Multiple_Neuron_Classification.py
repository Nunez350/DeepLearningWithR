import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import exp
from numpy.random import uniform

df = pd.read_csv("iris.csv")
df = df.drop([df.columns[0]], axis=1)
df = df.reset_index(drop=True)

# Multiple Neuron Classification Normalized by Softmax Function
def multiple_neuron_classification(input_cols, target_col, eta, iteration):    
    # create a list to store historical accuracy
    accuracy = []
    
    # declear inputs
    x = np.array(input_cols) # input dataset (dimension = number of individuals * number of features)
    X = np.insert(x, x.shape[1], 1, axis=1) # Matrix X: insert last column being the coefficients for bias (dimension = number of individuals * (number of features + 1))
    
    # build target values
    arr_target = np.array(target_col)
    ls = np.array(list(set(arr_target)))
    t = np.zeros([len(arr_target), len(ls)])
    for ix in range(len(ls)):
        t[np.where(arr_target == ls[ix]), ix] = 1
    # get index of ones by row for comparison with output y to compute accuracy of classification
    t_one_ix = np.argmax(t, axis=1)
    
    # initialize random weights (dimension of weights = features * neurons)
    weights = uniform(low=1e-3, high=1e-2, size=[x.shape[1], t.shape[1]])
    # initialize random bias (dimension of bias = 1 * neurons)
    bias = uniform(size=t.shape[1])
    # Matrix W: insert last row for bias (dimension = (number of weights + 1) * number of neurons)
    W = np.insert(weights, weights.shape[0], bias, axis=0)
    # create a collection to store historical calculated weights & append the initially generated weights
    weights_variation = {}
    for neuron_ix in range(weights.shape[1]):
        neuron_name = "Neuron_" + str(neuron_ix + 1)
        weights_variation[neuron_name] = {}
        for weight_ix in range(weights.shape[0]):
            weight_name = "W" + str(weight_ix + 1)
            weights_variation[neuron_name][weight_name] = []
            weights_variation[neuron_name][weight_name].append(weights[weight_ix][neuron_ix])
    
    # learning process
    for i in range(iteration):
        # neuron activity/output: compute inner product of Matrix X & W
        inner_product = np.dot(X, W)
        # compute y by normalizing inner_product of X & W based on softmax function
        y = exp(inner_product) / exp(inner_product).sum(axis=1).reshape(t.shape[0], 1)
        # compute accuracy of classification
        y_one_ix = np.argmax(y, axis=1)
        acc = np.where(np.equal(t_one_ix, y_one_ix) == True)[0].size / t.shape[0]
        accuracy.append(acc)
        # initialize an empty array to translate
        out_class = np.zeros(t.shape)
        
        # compute error e
        e = t - y
        
        # update weights & bias
        for neuron in range(W.shape[1]):
            neuron_name = "Neuron_" + str(neuron + 1)
            # update weights
            for weight in range(W.shape[0] - 1):
                weight_name = "W" + str(weight + 1)
                W[weight, neuron] = W[weight, neuron] - eta * ((-e[:, neuron] * X[:, weight]).sum() / e.shape[0])
                weights_variation[neuron_name][weight_name].append(W[weight, neuron])
            # update bias
            W[W.shape[0] - 1, neuron] = W[W.shape[0] - 1, neuron] - eta * ((-e[:, neuron]).sum() / e.shape[0])
    #return(weights_variation)
    return(accuracy, weights_variation)

def accuracy_variation(accuracy_array):
    curve, = plt.plot(accuracy_array, color='k', linewidth=1)
    plt.xlim([-10, 1000])
    plt.title("Accuracy Variation: Multiple Neuron Classification")
    plt.xlabel("Iteration")
    plt.ylabel("Accuracy")
    plt.show()

def weights_variation(weight_variation):
    line_style = ["-", "--", "-.", ":"]
    color_platte = ["r", "g", "b"]

    for neuron in weight_variation.keys():
        neuron_ix = int(neuron.replace("Neuron_", "")) - 1
        for weight in weight_variation[neuron].keys():
            weight_ix = int(weight.replace("W", "")) - 1
            curve, = plt.plot(weight_variation[neuron][weight],
                              label = neuron + ": " + weight,
                              color = color_platte[neuron_ix],
                              linestyle = line_style[weight_ix], linewidth=1)
            
    plt.legend(loc=2)
    plt.title("Variation of Weights: Multiple Neuron Classification")
    plt.xlabel("Iteration")
    plt.ylabel("Weight Value")
    plt.show()

if __name__ == "__main__":
    outputs = multiple_neuron_classification(df[df.columns[:4]], df[df.columns[4]], 0.1, 1000)
    weights_variation(outputs[1])
    accuracy_variation(outputs[0])
