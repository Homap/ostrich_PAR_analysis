from popgene.harmonic import harmonic_number
# **********************************
# Function for Watterson's theta
# **********************************
def theta_w(sn, n):
        """
        sn: number of segregating sites
        n: number of gene copies
        """
        end = n - 1
        har_num = harmonic_number(end)
        theta_watt = sn / har_num
        return theta_watt 