import numpy as np
from scipy.linalg import block_diag
from scipy.linalg import eigh
from scipy.sparse import dia_matrix
from scipy.sparse import spdiags  # Essentially not needed
import matplotlib.pyplot as plt
import math
import os


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def u_submatrix(site_matrix):
    # where n are eigenvalues and m are eigenvectors in a matrix form
    # n, m = LA.eigh(site_matrix)

    n, m = eigh(site_matrix)

    m = m[..., ::-1]
    m = m[::-1]

    return m, n


def u_matrix(ham_mol, ham_mollead, ham_lead, ham_em):  # where sub_U_mol_lead is extended molecule part

    u_lead, n_l = u_submatrix(ham_lead)
    print "lead eigenvals", n_l
    if split_type == "3split":
        u_em, n_em = u_submatrix(ham_em)
        print "em eigenvals", n_em
        uni_tran_mat = block_diag(u_lead, u_em, u_lead)
    elif split_type == "5split":
        u_m, n_m = u_submatrix(ham_mol)
        print "mol eigenvals", n_m
        u_ml, n_ml = u_submatrix(ham_mollead)
        print "ml eigenvals", n_ml
        uni_tran_mat = block_diag(u_lead, u_ml, u_m, u_ml, u_lead)
    else:
        print "u_matrix function has an error."

    uni_tran_mat_t = uni_tran_mat.transpose()
    uni_tran_mat_t = uni_tran_mat_t.conjugate()
    uni_tran_mat, uni_tran_mat_t = np.mat(uni_tran_mat), np.mat(uni_tran_mat_t)
    return uni_tran_mat, uni_tran_mat_t


def site_to_state_ham(site_ham, ham_mol, ham_mollead, ham_lead, ham_ext):
    u, u_t = u_matrix(ham_mol, ham_mollead, ham_lead, ham_ext)
    state_ham = np.dot(np.dot(u_t, site_ham), u)
    state_ham = np.array(state_ham)
    # print state_ham

    state_ham = np.mat(state_ham)
    return state_ham, u, u_t


# Molecule-lead couplings

# Needed for substitutions to be possible. I here choose specific couplings between the molecule itself and leads.
# as of now it is equal to beta

# IT GENERATES FULL HAMILTONIAN IN STATE REPRESENTATION. NEXT STEP: USE IT TO MAKE DENSITY MATRIX AND PROPAGATE IT.

# Split ham_state into ham_lead_left, ham_em and ham_lead_right. Then populate seperatly and combine

# number two in pop's  is because all molecular orbitals are doubly occupied.

# Perform reduction from full fermi function to if E-chem_pot < 0 then set 2.0 else 0.0


def step_function(energy, chem_pot):
    if abs(energy - chem_pot) < 10 **(-14):
        pop = 1.0
    elif energy - chem_pot < 0:
        pop = 2.0
    else:
        pop = 0.0
    return pop


def fermi_dirac_dist(energy, chem_pot, temp):
    k_b = 8.6173303 * 10 ** (-5)  # eV/K
    return 2.0 / (math.exp((energy - chem_pot) / (k_b * temp)) + 1.0)  # For low temperatures overflow error occurs.


def fermi_dirac_split(state_ham, temp, fermi_energy, bias_voltage, size_lead, size_wo_rl):
    energies = np.diag(state_ham, )
    populations = []

    chem_pot_ll = fermi_energy + 0.5 * bias_voltage
    chem_pot_rl = fermi_energy - 0.5 * bias_voltage

    for en in range(len(energies)):  # change to enumerate
        energy = float(energies[en])
        if en < size_lead:
            # Change of sign of chem_pot because of its symmetry around 0. E_F = 0 eV atm.
            if temp == 0:
                pop = step_function(energy, chem_pot_ll)
                populations.append(pop)
            else:
                pop = fermi_dirac_dist(energy, chem_pot_ll, temp)
                populations.append(pop)
        elif en >= size_wo_rl:
            if temp == 0:
                pop = step_function(energy, chem_pot_rl)
                populations.append(pop)
            else:
                pop = fermi_dirac_dist(energy, chem_pot_rl, temp)
                populations.append(pop)
        elif size_l <= en < (size_l + size_ml) and (split_type == "5split"):
            if temp == 0:
                pop = step_function(energy, chem_pot_ll)
                populations.append(pop)
            else:
                pop = fermi_dirac_dist(energy, chem_pot_ll, temp)
                populations.append(pop)
        elif (size_l + size_m + size_ml <= en < size_wo_rl) and (split_type == "5split"):
            if temp == 0:
                pop = step_function(energy, chem_pot_rl)
                populations.append(pop)
            else:
                pop = fermi_dirac_dist(energy, chem_pot_rl, temp)
                populations.append(pop)
        elif (size_l + size_ml) <= en < (size_l + size_ml + size_m) and (split_type == "5split"):
            if mol_loc_pot == "flat":
                pop = step_function(energy, 0.0)
                populations.append(pop)

            elif start_pot == "ramp":
                # print "LOL"
                loc_vol = - loc_vols[en - size_l - size_ml]
                if temp == 0:
                    pop = step_function(energy, loc_vol)
                    populations.append(pop)
                else:
                    pop = fermi_dirac_dist(energy, loc_vol, temp)
                    populations.append(pop)
            elif start_pot == "flat":
                pop = step_function(energy, 0.0)
                populations.append(pop)
            else:
                print "5split in fermi dirac is wrong"

        elif size_l <= en < size_wo_rl and (split_type == "3split"):
            if mol_loc_pot == "flat":
                pop = step_function(energy, 0.0)
                populations.append(pop)

            elif mol_loc_pot == "ramp":
                loc_vol = - loc_vols[en - size_l]
                if temp == 0:
                    pop = step_function(energy, loc_vol)
                    populations.append(pop)
                else:
                    pop = fermi_dirac_dist(energy, loc_vol, temp)
                    populations.append(pop)
            else:
                print "Poisson potential is not implemented for whole extended molecule. Results will most probably " \
                      "be wrong."
        else:
            print "THIS SHOULD NOT BE PRINTED. SOMETHING IS WRONG WITH FERMI_DIRAC_SPLIT FUNCTION. FIX IT!"

    data = np.array(populations)
    offsets = ([0])
    density_matrix = dia_matrix((data, offsets), shape=(size_tot, size_tot)).toarray()
    return density_matrix


# A single step of RK4 propagation


def rungekutta4(hamiltonian, density_matrix, ts, size_lead, target_density_matrix, damping, lead_zeros,
                em_l_zeros, l_em_zeros, em_zeros):
    density_t = dens_t_der(hamiltonian, density_matrix, size_lead, size_em, target_density_matrix, damping, lead_zeros,
                           em_l_zeros, l_em_zeros, em_zeros)
    k1 = density_t * ts
    dens_k = density_matrix + k1 * 0.5
    k2 = dens_t_der(hamiltonian, dens_k, size_lead, size_em, target_density_matrix, damping, lead_zeros, em_l_zeros,
                    l_em_zeros, em_zeros) * ts
    dens_k = density_matrix + k2 * 0.5
    k3 = dens_t_der(hamiltonian, dens_k, size_lead, size_em, target_density_matrix, damping, lead_zeros, em_l_zeros,
                    l_em_zeros, em_zeros) * ts
    dens_k = density_matrix + k3
    k4 = dens_t_der(hamiltonian, dens_k, size_lead, size_em, target_density_matrix, damping, lead_zeros, em_l_zeros,
                    l_em_zeros, em_zeros) * ts
    density_matrix = density_matrix + (k1 + 2 * (k2 + k3) + k4) / 6
    #  print (k1+2*(k2+k3)+k4)/6
    return density_matrix


# Function to calculate time derivative of density matrix: Liouville von Neumann equation

# TAKE CARE hamiltonian and density_matrix have to be in the same representation(both state or both site)

def dens_t_der(hamiltonian, density_matrix, size_lead, size_em, target_density_matrix, damping, lead_zeros, em_l_zeros,
               l_em_zeros, em_zeros):
    planck = 6.582119 * 10 ** (-1)  # eV * fs
    if damping != 0.0:
        damp_term = damping_term(size_lead, size_em, density_matrix, target_density_matrix, damping, lead_zeros,
                                 em_l_zeros, l_em_zeros, em_zeros)
    else:
        tot = 2 * size_lead + size_em
        damp_term = np.zeros((tot, tot))

    commutator = np.dot(hamiltonian, density_matrix) - np.dot(density_matrix, hamiltonian)

    commutator = np.mat(commutator)

    result = (-1.0j / planck) * commutator + damp_term  # / planck  # /planck REMOVED PLANCK because gamma is in fs^-1
    # unit therefore bc of E= planck * omega -> gamma = omega, gamma = E/planck
    return result


def bond_current(pnnplus):
    electron_charge = 1.602177 * 10 ** (-19)
    planck = 6.582119 * 10 ** (-1)
    return (2.0 * beta_m * electron_charge * 10 ** 18 / planck) * pnnplus.imag
    # 10 ** 18 because of fs calc. Should be times 10**18, removed for readability/floating point


def calc_current(ham_site, density_matrix, timestep, n_step, size_l, size_em, damping, u, u_t, frequency, temperature,
                 plot):
    site_pops = []
    state_pops = []

    currents = []  # []for n in range(size_em)]
    pnns = []
    p_diags = []
    sums_e_mol = []
    sums_e_em = []
    # em_currents : Currents at each site of extended molecule. Used to produce a heatmap of local currents at each site
    # vs time
    em_currents = [[] for n in range (n_step)]

    if mol_loc_pot == "poisson" and (split_type == "5split"):
        zeros_vector = np.zeros(size_m)
        zeros_vector[0] = 0.5 * bias_voltage
        zeros_vector[-1] = - 0.5 * bias_voltage
        zeros_vector *= epsilon / (atom_dist ** 2)
        beta_list = np.ones(size_m - 1)
        alfa_list = - 2 * np.ones(size_m)
        pot_matrix = - epsilon / (atom_dist ** 2) * tridiagonal(beta_list, alfa_list, beta_list)

    elif mol_loc_pot == "poisson" and (split_type == "3split"):
        zeros_vector = np.zeros(size_em)
        zeros_vector[0] = 0.5 * bias_voltage
        zeros_vector[-1] = - 0.5 * bias_voltage
        zeros_vector *= epsilon / (atom_dist ** 2)
        beta_list = np.ones(size_em - 1)
        alfa_list = - 2 * np.ones(size_em)
        pot_matrix = - epsilon / (atom_dist ** 2) * tridiagonal(beta_list, alfa_list, beta_list)
    else:
        print "Not poisson"

    ham_site_start = ham_site
    ham_state_start = np.dot(np.dot(u_t, ham_site_start), u)

    density_matrix = np.array(density_matrix)

    target_density_Ll = np.mat([row[0:size_l] for row in density_matrix[0:size_l]])
    target_density_Rl = np.mat([row[size_l + size_em:] for row in density_matrix[size_l + size_em:]])

    lead_zeros, em_l_zeros, l_em_zeros, em_zeros = all_zeros_matrices(size_l, size_em)
    matrix_target = damping * block_diag(target_density_Ll, em_zeros, target_density_Rl)

    l_to_pnn = size_l + size_em / 2
    # fs = 10**(-15)
    # eV*s.Change to eV*fs
    timestep = timestep  #
    # Coulomb per electron

    for ts in range(n_step):

        # Change to array so I can read elements of the matrix

        # Output elements of interest



        # Recover populations each 100'th step:
        if ts % 100 == 0:
            site_pops.append(np.diag(np.mat(density_matrix)))
            print "hey"
            state_dens = np.dot(np.dot(u_t, np.mat(density_matrix)), u)
            state_pops.append(np.diag(state_dens))
        density_matrix = np.array(density_matrix)
        sum_e_mol = np.diag(density_matrix)
        sum_e_em = np.diag(density_matrix)

        pnnplus = density_matrix[l_to_pnn - 1][l_to_pnn]

        p_dia = density_matrix[l_to_pnn - 1][l_to_pnn - 1]

        I = bond_current(density_matrix[l_to_pnn - 1][l_to_pnn])

        currents.append(I)
        pnns.append(pnnplus.imag)

        p_diags.append(p_dia.real)
        sums_e_mol.append(
            sum(sum_e_mol[size_l + size_ml:size_l + size_ml + size_m]).real)  # removed -1 from upper limit
        sums_e_em.append(sum(sum_e_em[size_l:size_l + size_em]).real)
        # Propagation using RK4
        global size_tot

        if heatmap == "yes":
            for site in range(size_em):
                em_site = size_l + site
                em_currents[ts].append(bond_current(density_matrix[em_site - 1][em_site]))
        elif heatmap != "no":
            print "Heatmap recievied incorrect user input."

        density_matrix = rungekutta4(ham_site, density_matrix, timestep, size_l, matrix_target, damping,
                                     lead_zeros,
                                     em_l_zeros, l_em_zeros, em_zeros)

        # Poisson magic
        if mol_loc_pot == "poisson" and (split_type == "5split"):
            rho_array = np.diag(density_matrix_site)[size_l + size_ml: size_l + size_ml + size_m]
            b_array = rho_array + zeros_vector

            x_vector = np.linalg.solve(pot_matrix, b_array)

            ham_site, loc_vols = poisson_iterate_exchange(ham_site, start_ham_m, x_vector)
        if mol_loc_pot == "poisson" and (split_type == "3split"):
            rho_array = np.diag(density_matrix_site)[size_l: size_l + size_em]
            b_array = rho_array + zeros_vector
            x_vector = np.linalg.solve(pot_matrix, b_array)

            ham_site, loc_vols = poisson_iterate_exchange(ham_site, start_ham_em, x_vector)
            # print "solving poisson eq!"
        # Some timedependent target matrix black magic:
        ts_sp = ts - starting_point

        if frequency != 0.0 and ts > starting_point:
            time_dep_pops = []
            bias_vol_time = chem_pot_time_dep(bias_voltage, frequency, ts_sp)  # minus starting point, so
            chem_pot_l = fermi_energy + 0.5 * bias_vol_time
            chem_pot_r = fermi_energy - 0.5 * bias_vol_time
            for ene in np.diag(ham_state_start)[0:size_l]:
                if temperature == 0.0:
                    pop = step_function(ene, chem_pot_l)
                    time_dep_pops.append(pop)
                else:
                    pop = fermi_dirac_dist(ene, chem_pot_l, temperature)
                    time_dep_pops.append(pop)
            for zero in np.diag(em_zeros):
                time_dep_pops.append(zero)
            for ene in np.diag(ham_state_start)[size_l + size_em:]:
                if temperature == 0.0:
                    pop = step_function(ene, chem_pot_r)
                    time_dep_pops.append(pop)
                else:
                    pop = fermi_dirac_dist(ene, chem_pot_r, temperature)
                    time_dep_pops.append(pop)
            tmat = dia_matrix((time_dep_pops, [0]), shape=(size_tot, size_tot)).toarray()
            matrix_target = damping * np.dot(np.dot(u, tmat), u_t)

        global gamma

        gamma = str(gamma)

        if plot == "plot" and (ts % 10 == 0):
            path = "Densities_vs_time" + type + gamma + split_type + "/Populations_SS_" + split_type + "_ts_" + str(ts)
            ensure_dir(path)
            real_pops = [item.real for item in np.diag(density_matrix)]
            filename = path
            plt.plot(real_pops)
            plt.ylabel("Population")
            plt.xlabel("Site number")
            plt.title("SS populations across the system")
            plt.savefig(filename)
            plt.close()

        if ts % 10 == 0:
            print ts

    return site_pops, state_pops, em_currents, currents, pnns, p_diags, sums_e_mol, sums_e_em


# Make Gamma matrix out of target density matrix.
# I split gamma matrix into three parts. leads and lead-lead coherences part, extended molecule-lead part and target density matrix part.

def damping_term(size_l, size_em, density_matrix, target_density_matrix, damping, lead_zeros, em_l_zeros, l_em_zeros,
                 em_zeros):
    # planck = 6.582119 * 10 ** (-1)  # Change to eV * fs
    density_matrix = np.array(density_matrix)

    # Lead-lead and lead-lead coherences. Capital letters for Left/Right
    Ll_density_matrix = np.mat([row[0:size_l] for row in density_matrix[0:size_l]])
    Rl_density_matrix = np.mat([row[size_l + size_em:] for row in density_matrix[size_l + size_em:]])
    # Right-Left lead and Left-Right lead coherences:
    RL_coh_density_matrix = np.mat([row[0:size_l] for row in density_matrix[size_l + size_em:]])
    LR_coh_density_matrix = np.mat([row[size_l + size_em:] for row in density_matrix[0:size_l]])

    lead_lead_coh = - damping * np.bmat(
        [[Ll_density_matrix, l_em_zeros, LR_coh_density_matrix], [em_l_zeros, em_zeros, em_l_zeros],
         [RL_coh_density_matrix, l_em_zeros, Rl_density_matrix]])

    # Extended molecule-lead coherences
    em_Ll_coh = np.mat([row[0:size_l] for row in density_matrix[size_l:size_l + size_em]])
    Ll_em_coh = np.mat([row[size_l:size_l + size_em] for row in density_matrix[0:size_l]])
    Rl_em_coh = np.mat([row[size_l:size_l + size_em] for row in density_matrix[size_l + size_em:]])
    em_Rl_coh = np.mat([row[size_l + size_em:] for row in density_matrix[size_l:size_l + size_em]])

    ext_mol_lead_coh = - 0.5 * damping * np.bmat(
        [[lead_zeros, Ll_em_coh, lead_zeros], [em_Ll_coh, em_zeros, em_Rl_coh], [lead_zeros, Rl_em_coh, lead_zeros]])
    result_mat = lead_lead_coh + ext_mol_lead_coh + target_density_matrix
    # print "leadleadcoh", lead_lead_coh
    # print "extmolLeadcoh", ext_mol_lead_coh
    # print "target dens mat", target_density_matrix
    return np.mat(result_mat)


def chem_pot_time_dep(chem_pot, frequency, time):
    return chem_pot * math.cos(
        frequency * 2.0 * math.pi * time)  # CARE WITH UNITS HERE. As theory -> time in fs , frequency in 1/ps


def all_zeros_matrices(size_l, size_em):
    lead_zeros = np.zeros((size_l, size_l))
    em_l_zeros = np.zeros((size_em, size_l))
    l_em_zeros = np.zeros((size_l, size_em))
    em_zeros = np.zeros((size_em, size_em))
    return lead_zeros, em_l_zeros, l_em_zeros, em_zeros


def ramp_func(v_bias, ham_mol):
    """

    :rtype: array
    """
    len_m = len(np.diag(ham_mol))
    loc_vols = []
    for i, alfa in enumerate(np.diag(ham_mol)):
        loc_vol = (((i + 1) * v_bias) / (len_m + 1)) + 0.5 * v_bias
        alfa_td = alfa - loc_vol
        ham_mol[i][i] = alfa_td
        loc_vols.append(loc_vol)
    return ham_mol, loc_vols

def benzene_ham():
    ham_mol = np.array(ham_m)
    if benzene_type == "para":
        ham_mol[2][3] = 0.0
        ham_mol[3][2] = 0.0
        ham_mol[0][3] = beta_m
        ham_mol[3][0] = beta_m
        ham_mol[2][5] = beta_m
        ham_mol[5][2] = beta_m
    elif benzene_type == "meta":
        ham_mol[4][3] = 0.0
        ham_mol[3][4] = 0.0
        ham_mol[0][4] = beta_m
        ham_mol[4][0] = beta_m
        ham_mol[0][5] = beta_m
        ham_mol[5][0] = beta_m
    elif benzene_type == "ortho":
        ham_mol[0][5] = beta_m
        ham_mol[5][0] = beta_m
    else:
        print "Wrong type specification"

    return ham_mol

def tridiagonal(beta_off_minus, alfa_a, beta_off_plus):
    return np.diag(beta_off_minus, -1) + np.diag(alfa_a, 0) + np.diag(beta_off_plus, 1)


def poisson_iterate_exchange(ham_site, start_ham_molecule, solution_vector):
    new_loc_vols = []
    ham_site, start_ham_molecule = np.array(ham_site), np.array(start_ham_molecule)
    if split_type == "5split":
        for i in range(size_m):
            ham_site[i + size_l + size_ml][i + size_l + size_ml] = start_ham_molecule[i][i] - solution_vector[i]
            new_loc_vols.append(solution_vector[i])

    elif split_type == "3split":
        for i in range(size_em):
            ham_site[i + size_l][i + size_l + size_em] = start_ham_molecule[i][i] - solution_vector[i]
            new_loc_vols.append(solution_vector[i])
    else:
        print "error in poisson iterate"
    return ham_site, new_loc_vols

# SIZES
size_m = int(raw_input("Size of molecule matrix: "))
size_ml = int(raw_input("Size of molecule-lead matrices: "))
size_l = int(raw_input("Size of electrode matrices: "))

# SYSTEM DESCRIPTORS
benzene = raw_input("Is molecule benzene? yes or no: ")
if benzene == "yes":
    benzene_type = raw_input("Insert whether it is para, meta or ortho bonded: ")  # Write "para", "meta" or "ortho"
elif benzene != "no":
    print "Please write 'no'"

direct_coupling = float(raw_input("Size of direct molecule-lead coupling: "))  # eV
coupling = float(raw_input("Insert Molecule-Lead coupling V: "))
beta_l = float(raw_input("Insert hopping matrix element in leads: "))
beta_ml = beta_l
beta_m = float(raw_input("Insert hopping matrix element in molecule: "))
bias_voltage = float(raw_input("Insert bias voltage in V: "))
fermi_energy = float(raw_input("Insert fermi energy in eV: "))
temperature = float(raw_input("Insert electronic temperature: "))
mol_loc_pot = raw_input("flat, ramp or poisson potential: ")
split_type = raw_input("System splitting in 3 or 5 parts. Write 3split or 5split: ")
if mol_loc_pot == "poisson":
    epsilon = float(raw_input("Molecular permativity: "))
    atom_dist = float(raw_input("Interatomic distance in atomic units: ")) * 5.29177 * 10 ** (-11)
    start_pot = raw_input("Start Potential type flat or ramp: ")
elif mol_loc_pot == "ramp":
    start_pot = "ramp"
else:
    start_pot = "zero"



# SIMULATION DESCRIPTORS
number_steps = int(raw_input("Number of steps: "))
time = range(number_steps)
gamma = float(raw_input("Size of driving factor: "))  # SIZE OF DRIVING FACTOR
frequency = float(raw_input("Insert non-zero frequency in fs^-1 unit if AC current is desired, else insert 0.0: "))
starting_point = 0
if frequency > 0.0:
    starting_point = int(raw_input("Write after how many timesteps does frequency start: "))
type = raw_input("Calculation type: ")  # Characterstic for a series of calculations
plot = raw_input("plot or test: ")
heatmap = raw_input("Produce heatmap across extended molecule. yes or no: ")


# PROMPT FOR TYPE IN MOLECULE. RAMP / POISSON / FLAT

size_em = size_m + size_ml * 2

size_tot = size_em + 2 * size_l
size_wo_rl = size_l + size_em

diags = np.array([-1, 1])  # Offdiagonal hopping matrix elements

diag_long_up = np.array([size_l - 1])
diag_long_down = np.array([-size_l + 1])
diag_short_up = np.array([size_m + size_ml * 2 - 1])
diag_short_down = np.array([-size_m - size_ml * 2 + 1])

# HAM, with possibility to be alter, leads, molecule and 9extended molecule part. Change to simpler: Change one element
# at a time
data_m = np.array([[beta_m for s in range(size_m)], [beta_m for x in range(size_m)]])  # TEST CASE  MOLECULE BETA -> 1
data_ml = np.array([[beta_ml for s in range(size_ml)], [beta_ml for x in range(size_ml)]])
data_l = np.array([[beta_l for s in range(size_l)], [beta_l for x in range(size_l)]])
data_v = np.array([[coupling for s in range(size_l)]])

ham_m = spdiags(data_m, diags, size_m, size_m).toarray()
ham_ml = spdiags(data_ml, diags, size_ml, size_ml).toarray()
ham_l = spdiags(data_l, diags, size_l, size_l).toarray()

start_ham_m = ham_m

if mol_loc_pot != "flat" and (start_pot == "ramp") and (benzene != "yes"):
    if split_type == "5split":

        ham_m, loc_vols = ramp_func(bias_voltage, ham_m)
    else:
        print "Going to 3split"

if benzene == "yes":
    ham_m = benzene_ham()
else:
    print "proceeding with molecular wire"



# Combine ml and m elemenets into extenede molecule
# CORRECTION input direct couplings between molecule and lead.
ham_em = block_diag(ham_ml, ham_m, ham_ml)
# Try np.bmat instead ofblock_diag. FAIL
# ham_em = np.bmat(np.mat(ham_ml), np.mat(ham_m), np.mat(ham_ml))

ham_em = np.array(ham_em)
# print ham_em
# Direct couplings between molecule and lead

# Left lead to molecule
ham_em[size_ml][size_ml - 1] = beta_ml
ham_em[size_ml - 1][size_ml] = beta_ml

# Right lead to molecule
ham_em[size_ml + size_m][size_ml + size_m - 1] = coupling
ham_em[size_ml + size_m - 1][size_ml + size_m] = coupling

# Creation of 4 V non-squares matrices that couple extended molecule to leads.

V_ml = spdiags(data_v, diag_long_up, size_em, size_l).toarray()
V_lm = spdiags(data_v, diag_long_down, size_l, size_em).toarray()
V_mr = spdiags(data_v, diag_short_down, size_em, size_l).toarray()
V_rm = spdiags(data_v, diag_short_up, size_l, size_em).toarray()

# Also two lead_size times lead_size all zero matrices REMOVE LATER. THERE IS A ZEROS FUNCTION

lead_size_zeros = np.zeros((size_l, size_l))

# Hamiltonian in site representation:

if mol_loc_pot != "flat":
    if split_type == "3split":
        start_ham_em = ham_em
        ham_em, loc_vols = ramp_func(bias_voltage, ham_em)
    elif split_type == "5split":
        print "Proceeding with 5 split"
    else:
        print "Wrong split specification"

ham_l = np.mat(ham_l)
ham_em = np.mat(ham_em)

ham_site = np.bmat([[ham_l, V_lm, lead_size_zeros], [V_ml, ham_em, V_mr], [lead_size_zeros, V_rm, ham_l]])

ham_state, u, u_t = site_to_state_ham(ham_site, ham_m, ham_ml, ham_l, ham_em)

density_matrix_state = fermi_dirac_split(ham_state, temperature, fermi_energy, bias_voltage, size_l, size_wo_rl)

# Change to site representation

density_matrix_site = np.dot(np.dot(u, density_matrix_state), u_t)
density_matrix_site = np.array(density_matrix_site)

# START

# EXECUTE

site_populations, state_populations, em_currents, currents, pnns, p_diag, sums_e_mol, sums_e_em = calc_current(ham_site,
                                                                                                    density_matrix_site,
                                                                                                    1.0,
                                                                                                    number_steps,
                                                                                                    size_l,
                                                                                                    size_em,
                                                                                                    gamma, u,
                                                                                                    u_t,
                                                                                                    frequency,
                                                                                                    temperature,
                                                                                                    plot)

if heatmap == "yes":
    heatimage = plt.imshow(em_currents, cmap='hot', interpolation='nearest')
    cbar = plt.colorbar(heatimage)
    plt.ylabel("Time(fs)")
    plt.xlabel("Extended Molecule site no.")
    cbar.set_label("Currents (mA)")
    if plot == "test":
        plt.show()
        plt.close()
    elif plot == "plot":
        heatname = type + ".png"
        plt.savefig(heatname)
        plt.close()

if plot == "plot":

    if heatmap == "yes":
        with open("Currents_vs_em_site"+type+"_.txt", 'w') as em:
            for t in em_currents:
                for cur in t:
                    cur = str(cur)
                    em.write(cur)
                    em.write(" ")
                em.write("\n")

    with open("Energies_state_ham_diag" + type + "_.txt", 'w') as f:
        for energy in np.diag(np.array(ham_state)):
            energy = str(energy)
            f.write(energy)
            f.write(" ")

    with open("Populations_site" + type + "_.txt", 'w') as g:
        for populations in site_populations:
            for pop in populations:
                pop = str(pop.real)
                g.write(pop)
                g.write(" ")
            g.write("\n")

    with open("Populations_state" + type + "_.txt", 'w') as h:
        for populations in state_populations:
            for pop in populations:
                pop = str(pop.real)
                h.write(pop)
                h.write(" ")
            h.write("\n")

    with open("Currents_1fs" + type + str(gamma) + ".txt", "w") as f:  # Currents_pnns_pdiag_mol_em_500_01K.txt"
        f.write("Current" + " " + "pnn" + " " + "p_diag" + " " + "sum_e_mol" + " " + "sum_e_em" + "\n")
        for current, pnn, p_diag, sum_e_mol, sum_e_em in zip(currents, pnns, p_diag, sums_e_mol, sums_e_em):
            current = str(current)
            pnn = str(pnn)
            p_diag = str(p_diag)
            sum_e_mol = str(sum_e_mol)
            sum_e_em = str(sum_e_em)
            f.write(current + " " + pnn + " " + p_diag + " " + sum_e_mol + " " + sum_e_em)
            f.write("\n")

    plt.plot(time, currents)
    plt.xlabel("Time (fs)")
    plt.ylabel("Current (mA)")
    plt.title("Damping " + str(gamma))
    plt.savefig("Currents_1fs_" + str(gamma) + type + ".png")
    plt.close()

    plt.plot(time, sums_e_mol)
    plt.xlabel("Time (fs)")
    plt.ylabel("Number of e- at molecule")
    plt.title("Damping " + str(gamma))
    plt.savefig("sums_e_mol_" + str(gamma) + "_1fs" + type + ".png")
    plt.close()

    plt.plot(time, sums_e_em)
    plt.xlabel("Time (fs)")
    plt.ylabel("Number of e- at extended molecule")
    plt.title("Damping " + str(gamma))
    plt.savefig("sums_e_em_" + str(gamma) + "_1fs_" + type + ".png")
    plt.close()

    plt.plot(time, pnns)
    plt.xlabel("Time (fs)")
    plt.ylabel("Coherence")
    plt.title("Damping " + str(gamma))
    plt.savefig("coherence_" + str(gamma) + "_1fs_" + type + ".png")
    plt.close()
elif plot == "test":
    plt.plot(currents)
    plt.show()
