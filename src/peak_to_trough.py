"""
Peak to trough analysis
"""

import matplotlib.pyplot as plt
import random

random.seed(321)

def drawdown(prices):
    prevmaxi = 0
    prevmini = 0
    maxi= 0

    for i in range(len(prices))[1:]:
        if prices[i] >= prices[maxi]:
            maxi = i
        else:
            if (prices[maxi] - prices[i]) > (prices[prevmaxi] - prices[prevmini]):
                prevmaxi = maxi
                prevmini = i
    return [[prevmaxi, prevmini], [prices[prevmaxi], prices[prevmini]]]

def drawdown_reverse(prices):
    prevmaxi = 0
    prevmini = 0
    mini=0

    for i in range(len(prices))[1:]:
        if prices[i] <= prices[mini]:
            mini = i            
        else:
            if (prices[i] - prices[mini]) > (prices[prevmaxi] - prices[prevmini]):
                prevmaxi = i
                prevmini = mini
    return [[prevmaxi, prevmini], [prices[prevmaxi], prices[prevmini]]]

def time_series_data(alpha, beta, steps, sigma):
    y = [0]
    l = [0]
    b = [0]

    #err = random.gauss(0, sigma)

    for i in range(steps):
        y.append(l[i] + b[i] + random.gauss(0, sigma))
        l.append(l[i] + b[i] + alpha * random.gauss(0, sigma))
        b.append(b[i] + alpha * beta * random.gauss(0, sigma))

    return y

length_series = 100
ts_data = time_series_data(1.5, 1, length_series, 1.5)
peak_to_trough = drawdown(ts_data)
trough_to_peak = drawdown_reverse(ts_data[peak_to_trough[0][1]:])
#trough_to_peak = drawdown_reverse(ts_data)
print trough_to_peak

plt.plot(range(length_series+1), ts_data)
plt.plot(peak_to_trough[0], peak_to_trough[1])
plt.plot(trough_to_peak[0], trough_to_peak[1])
plt.show()
