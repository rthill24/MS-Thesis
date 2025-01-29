y = 0
tolerance = 0.0001
while True:
    y_new = (4*(y**2)+6)/10
    if abs(y_new-y) < tolerance:
        break
    y = y_new
print(y_new)