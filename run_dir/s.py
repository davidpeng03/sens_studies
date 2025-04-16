
import matplotlib.pyplot as plt
import datetime

now = datetime.now()
plotSuffix = now.strftime("%Y%m%d_%H_%M_S")
plt.plot([1, 2, 3], [4, 5, 6])
plt.show()

print("...")
a = input()