# Configuration template
ConfigName = "Configuration I"

# In this file, you can update the parameters that define the geometry and the materials of the lug we're designing.
# Please make a copy of this template when adjusting and give it a new name. You can test your configuration
# by changing the import statement in the Master file.

# Spacecraft skin:
t3      = 2     * 10**(-3)      #[m]        Thickness
Eskin   = 72    * 10**(9)       #[Pa]       Young Modulus
Gskin   = 27    * 10**(9)       #[Pa]       Shear Modulus

# Back Plate:
t2      = 5     * 10**(-3)      #[m]        Thickness
L       = 250   * 10**(-3)      #[m]        Lenght of the back plate
w       = 100   * 10**(-3)      #[m]        Height of the back plate
D2      = 4     * 10**(-3)      #[m]        Diameter of bolts/rivets
n       = 4                     #           Number of bolts in the back plate (total)
e1      = 5     * 10**(-3)      #[m]        distance edge to center of a hole on the y-axis (along w dimension)
e2      = 5     * 10**(-3)      #[m]        distance edge to center of a hole on the x-axis (along L dimension)

# Flange
s1      = 100   * 10**(-3)      #[m]        Length of the flange
s2      = 70    * 10**(-3)      #[m]        Center position of the pin hole
D1      = 30    * 10**(-3)      #[m]        Diameter of the pin hole
h       = 50    * 10**(-3)      #[m]        Inner dimension between two flanges
N       = 2                     #           Number of flanges

# Material lug
E       = 72    * 10**(9)       #[Pa]       Young Modulus
G       = 27    * 10**(9)       #[Pa]       Shear Modulus