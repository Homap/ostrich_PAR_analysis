# **********************************
# Function for harmonic number
# **********************************
def harmonic_number(n):
        s_har = 0
        for i in range(1, n+1):
                s_har = s_har + 1/i
        return s_har