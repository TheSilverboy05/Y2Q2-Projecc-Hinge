# Configuration template
ConfigName = "Configuration I"

# In this file, you can update the parameters that define the geometry and the materials of the lug we're designing.
# Please make a copy of this template when adjusting and give it a new name. You can test your configuration
# by changing the import statement in the Master file.

# Main forces on each lug
Fx      = 100   * 10**(0)       #[N]        Main load in x-direction
Fy      = 100   * 10**(0)       #[N]        Main load in y-direction
Fz      = 100   * 10**(0)       #[N]        Main load in z-direction

# Main Moments on each lug
Mx      = 0     * 10**(0)       #[Nm]       Main moment in x-direction
My      = 2     * 10**(0)       #[Nm]       Main moment in y-direction
Mz      = 2     * 10**(0)       #[Nm]       Main moment in z-direction

# Spacecraft skin:
t3      = 2     * 10**(-3)      #[m]        Thickness
Eskin   = 72    * 10**(9)       #[Pa]       Young Modulus
Gskin   = 27    * 10**(9)       #[Pa]       Shear Modulus

# Lug Distance & Solar panel dimension
H       = 500   * 10**(-3)      #[m]        Distance between the lugs
W       = 1200  * 10**(-3)      #[m]        Length of the deployed solar panel

# Back Plate:
t2      = 5     * 10**(-3)      #[m]        Thickness
L       = 250   * 10**(-3)      #[m]        Lenght of the back plate
w       = 75    * 10**(-3)      #[m]        Height of the back plate
D2      = 4     * 10**(-3)      #[m]        Diameter of bolts/rivets
n       = 4                     #           Number of bolts in the back plate (total)
e1      = 5     * 10**(-3)      #[m]        distance edge to center of a hole on the y-axis (along w dimension)
e2      = 5     * 10**(-3)      #[m]        distance edge to center of a hole on the x-axis (along L dimension)

# Flange
t1      = 3     * 10**(-3)      #[m]        Thickness of the flange
s1      = 100   * 10**(-3)      #[m]        Length of the flange
s2      = 70    * 10**(-3)      #[m]        Center position of the pin hole
D1      = 70    * 10**(-3)      #[m]        Diameter of the pin hole
h       = 50    * 10**(-3)      #[m]        Inner dimension between two flanges
N       = 2                     #           Number of flanges

# Material lug
E       = 72    * 10**(9)       #[Pa]       Young Modulus
G       = 27    * 10**(9)       #[Pa]       Shear Modulus


#Materials Bolts


class Material:
    def __init__(self, id, name, tensile_strength, shear_strength, density):

        self.id = id 
        self.name = name 
        self.tensile_strength = tensile_strength
        self.shear_strength = shear_strength
        self.density = density

    def __repr__(self):
        """
        Provide a string representation of the material for easy viewing.
        """
        return (f"Material {self.id}: {self.name}\n"
                f"  Tensile Yield Strength: {self.tensile_strength} MPa\n"
                f"  Shear Yield Strength: {self.shear_strength} MPa\n"
                f"  Density: {self.density} kg/m³\n")


# Defining the 4 materials
materials = [
    Material(0, "7075 T6 Alluminium Alloy", 483*10**6, 331*10**6, 2810),   # properties
    Material(1, "2014 T6 Alluminium Alloy", 400*10**6, 290*10**6, 2800),
    Material(2, "SAE-AISI 4340 Steel", 470*10**6, 430*10**6, 7800),
    Material(3, "Aged Grade 250 Maraging Steel", 1740*10**6, 1060*10**6, 8200),
]

# Display the materials
for material in materials:
    print(material)
    
print(materials[0].shear_strength) #test test use this if you need any feature
