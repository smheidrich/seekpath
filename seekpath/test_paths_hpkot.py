from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import str
from builtins import range
import unittest
from .util import atoms_num_dict


def simple_read_poscar(fname):
    with open(fname) as f:
        lines = [l.partition('!')[0] for l in f.readlines()]

    alat = float(lines[1])
    v1 = [float(_) * alat for _ in lines[2].split()]
    v2 = [float(_) * alat for _ in lines[3].split()]
    v3 = [float(_) * alat for _ in lines[4].split()]
    cell = [v1, v2, v3]

    species = lines[5].split()
    num_atoms = [int(_) for _ in lines[6].split()]

    next_line = lines[7]
    if next_line.strip().lower() != 'direct':
        raise ValueError(
            "This simple routine can only deal with 'direct' POSCARs")
    # Note: to support also cartesian, remember to multiply the coordinates
    # by alat

    positions = []
    atomic_numbers = []
    cnt = 8
    for el, num in zip(species, num_atoms):
        atom_num = atoms_num_dict[el.capitalize()]
        for at_idx in range(num):
            atomic_numbers.append(atom_num)
            positions.append([float(_) for _ in lines[cnt].split()])
            cnt += 1

    return (cell, positions, atomic_numbers)

class TestTransformationMatrix(unittest.TestCase):
    """
    Test the calculation of the transformation matrix
    """

    # List of transformations to test on various systems
    transformations_to_test = [
        # Start checking the identity
        [[ 1.,  0.,  0.], 
         [ 0.,  1.,  0.],
         [ 0.,  0.,  1.]],
        # Swap two axes a-c (this also changes the handedness!)
        [[ 0.,  0.,  1.], 
         [ 0.,  1.,  0.],
         [ 1.,  0.,  0.]],
        # Swap two axes a-b
        [[ 0.,  1.,  0.], 
         [ 1.,  0.,  0.],
         [ 0.,  0.,  1.]],
        # Swap two axes b-c
        [[ 1.,  0.,  0.], 
         [ 0.,  0.,  1.],
         [ 0.,  1.,  0.]],
        # I only swap 1 axis
        [[ -1.,  0.,  0.], 
          [ 0.,  1.,  0.],
          [ 0.,  0.,  1.]],
        # I only swap the second axis
        [[ 1.,  0.,  0.], 
          [ 0.,  -1.,  0.],
          [ 0.,  0.,  1.]],
        # I only swap the third axis
        [[ 1.,  0.,  0.], 
          [ 0.,  1.,  0.],
          [ 0.,  0.,  -1.]],
        # I swap all three axis
        [[ -1.,  0.,  0.], 
          [ 0.,  -1.,  0.],
          [ 0.,  0.,  -1.]],
        # I also swap 1 axis
        [[ -1.,  0.,  0.], 
          [ 0.,  0.,  -1.],
          [ 0.,  -1.,  0.]],
        # I also swap 1 axis, but with a different sign convention
        [[ -1.,  0.,  0.], 
          [ 0.,  0.,  1.],
          [ 0.,  1.,  0.]],
        ]

    def _is_same_transf_matrix(self, m1, m2, conv_lattice, thr=1.e-8):
        """
        Check if the two matrices m1 and m2 (that must have the same shape)
        represent the same transformation.
        This consists in two steps:

        1. check if the two matrices are almost equal (in the sense defined
           below in point 2) for at least one change in sign of any pair of 
           axes that would preserve vector lengths, angles and right-handedness
        2. the two matrices (with possibly a change in sign of two axes as
           explained in point 1) are almost equal, within the given 
           threshold ``thr`` (i.e., the the mean of the absolute differences
           of each component is smaller than ``thr``)

        Explanation for point 1: a swap of sign of two vectors ``v_1`` and 
        ``v_2``is equivalent to a rotation of 180 degrees around the 
        orthogonal axis, if the third vector ``v_3`` is orthogonal to both
        ``v_1`` and ``v_2``. In this case, the length of the three vectors 
        and their mutual angles are unchanged, as well as their handedness.
        In this case, this rotation can be arbitrary put in the transformation 
        matrix, or in the rotation of the Cartesian axes 
        (presenving the right-handedness of the three vectors).

        The separation of the global transformation in a transformation matrix
        mixing the vectors, and a rotation of the Cartesian coordinates, is
        therefore partially arbitrary (in these cases with at leaast one 
        vector orthogonal to the other two), and we need that at least one such
        matrix is equal.

        We therefore apply the transformation, check that the six lattice
        parameters did not change, and in this case consider also this matrix.
        For this reason, we also need the standardized lattice ``conv_lattice``.
        Pass, there, the standardized (conventional) lattice, 
        because it is also refined.
        """
        import numpy as np
        m1arr = np.array(m1, dtype=float)
        m2arr = np.array(m2, dtype=float)
        lattice = np.array(conv_lattice, dtype=float)

        # All possible sign swaps (except the identity, always added belop)
        sign_swap_matrices = [
            np.array([
                [1.,0.,0.],
                [0.,-1.,0.],
                [0.,0.,-1.]]),
            np.array([
                [-1.,0.,0.],
                [0.,1.,0.],
                [0.,0.,-1.]]),
            np.array([
                [-1.,0.,0.],
                [0.,-1.,0.],
                [0.,0.,1.]])        
            ]

        # Identity must always be checked
        valid_swaps = [
            np.array([
            [1.,0.,0.],
            [0.,1.,0.],
            [0.,0.,1.]])
            ]

        for sign_swap in sign_swap_matrices:
            lattice_with_swapped_vectors = np.matmul(sign_swap, lattice)

            # cell_matrix * transpose(cell_matrix): on the diagonal: 
            # square of vector lengths; out of diagonal: it's related to the 
            # cosine of the angles
            # We want that cell parameters do not change, equivalent to 
            # checking that these products do not change
            lengths_and_angles_orig = np.matmul(
                lattice, np.array(lattice).T)
            lengths_and_angles_swapped = np.matmul(
                lattice_with_swapped_vectors, lattice_with_swapped_vectors.T)

            # Add to the valid swaps only if the vector swap is equivalent
            # to a rotation, i.e. doesn't change lengths and vectors.
            if np.mean(np.abs(
                lengths_and_angles_swapped - lengths_and_angles_orig)) < thr:
                valid_swaps.append(sign_swap)

        # Check that at least one valid swap produces the same transformation
        # matrix
        # The swap is multiplied on the right: in fact, 
        # we just checked above that the std_lattice is equivalent also if 
        # replaced with ``swap * std_lattice``. Since we also have the relation
        # input_lattice = T * std_latt * R, 
        # we have also
        # input_lattice = T * swap * std_latt * R2, 
        # (with R2 a rotation counteracting the axis swap)
        # so the new transformation matrix is (T * swap)
        equal = False
        for valid_swap in valid_swaps:
            m2arr_with_swap = np.matmul(
                np.array(m2arr), np.array(valid_swap))
            mean_difference = np.mean(np.abs(m1arr - 
                m2arr_with_swap))
            equal = equal or (mean_difference < thr)

        self.assertTrue(equal, 
            "The two matrices are different (for any valid swap of signs of "
            "pairs of axes): {} vs {}".format(
                m1arr, m2arr
                )
            )

    def _test_vector_permutations(self, system, basic_transf):
        """
        Test all permutations defined in self.transformations_to_test,
        for the ``system`` given in input, checking that indeed we
        get the composition of ``basic_transf`` with the applied 
        transformation.

        :param system: must be a tuple in the form 
            (lattice, positions, atomic_numbers)
        :param basic_transf: must be a 3x3 array, containing the expected
            transformation for the given system (it is explicitly given as
            a further check - anyway the identity transformation is also 
            checked)
        """
        from seekpath import hpkot
        import numpy as np

        # This is already in the standard order a<b<c
        cell, positions, atomic_numbers = system

        for my_transf in self.transformations_to_test:

            assert abs(np.linalg.det(my_transf)) == 1, (
                "Only transformations that do not generate a supercell can "
                "be tested with this method, change the matrix in "
                "self.transformations_to_test (reason: it is more complex, "
                "otherwise, to apply the transformation to the positions)")

            tr_cell = np.matmul(my_transf, cell)
            # I must also apply the transformation to the atomic positions
            tr_positions = np.matmul(my_transf, np.array(positions).T).T

            system = (tr_cell, tr_positions, atomic_numbers)
            res = hpkot.get_path(system)
            conv_transformation_matrix = res['conv_transformation_matrix']
            # relation: initial_cell = T[basic_transf] * std_latt * R"
            # if I multiply my_transf by initial_cell (to get tr_cell)
            # I have: 
            # my_transf*initial_cell = [my_transf*basic_transf] * std_latt * R
            # and my_transf*basic_transf should therefore be equal to the 
            # transformation matrix I get from the code
            if res['conv_transformation_matrix'] is None:
                raise AssertionError(
                    "The conventional transformation matrix from "
                    "seekpath is none...")
            try:
                self._is_same_transf_matrix(
                    np.matmul(my_transf, basic_transf), 
                    conv_transformation_matrix, 
                    conv_lattice=res['conv_lattice'])
            except AssertionError as e:
                raise AssertionError(
                    "{} (my transformation matrix is {})".format(
                        e.message, my_transf))


    def test_ortho_swapped_axes(self):
        """
        Simple orthogonal case, axes swapped
        """
        from seekpath import hpkot

        cell = ((1,0,0),(0,0,2), (0,1.5,0))
        positions = ((0,0,0),)
        atomic_numbers = [1,]

        system = (cell, positions, atomic_numbers)
        res = hpkot.get_path(system, with_time_reversal=False)

        conv_transformation_matrix = res['conv_transformation_matrix']

        self._is_same_transf_matrix(conv_transformation_matrix, 
            [[ 1.,  0.,  0.],
             [ 0.,  0.,  1.],
             [ 0.,  1.,  0.]],
             conv_lattice=res['conv_lattice'])

    def test_permutations_ortho(self):
        """
        Simple orthogonal case,  do various permutations and check
        the conv_transformation_matrix
        """
        import numpy as np 
        cell = np.array([(1.,0.,0.),(0.,1.5,0.), (0,0.,2.)])
        positions = ((0,0,0),)
        atomic_numbers = [1,]

        # This is the transformation of the cell above
        # (identity, it should be already standardized
        basic_transf = np.diag([1,1,1])

        self._test_vector_permutations(
            system = (cell, positions, atomic_numbers),
            basic_transf=basic_transf
            )

    def test_permutations_supercell_ortho(self):
        """
        Simple orthogonal case, check permutations, but considering
        a 2x1x1 supercell
        """
        import numpy as np         
        # This is already in the standard order a<b<c
        cell = np.array([(1.,0.,0.),(0.,1.5,0.), (0,0.,2.)])
        positions = ((0,0,0),(0.5,0,0))
        atomic_numbers = [1,1]

        # This is the transformation of the cell above
        # (there is a 2 because it's a supercell)
        basic_transf = np.diag([2,1,1])

        self._test_vector_permutations(
            system = (cell, positions, atomic_numbers),
            basic_transf=basic_transf
            )

    def test_permutations_monoclinic(self):
        """
        Simple monoclinic case, do various permutations and check
        the conv_transformation_matrix.
        """
        import numpy as np         
        
        cell = np.array([(1.,0.,0.),(0.,2.2,0.),(-0.4,0.,1.5)])

        positions = ((0.,0.,0.),)
        atomic_numbers = [1]

        # This is the transformation of the cell above
        # (it's already standardized, so it is the identity)
        basic_transf = np.diag([1,1,1])

        self._test_vector_permutations(
            system = (cell, positions, atomic_numbers),
            basic_transf=basic_transf
            )

    def test_permutations_triclinic(self):
        """
        Simple orthogonal case, do various permutations and check
        the conv_transformation_matrix (it should be the inverse, as I
        start from a cell that is already standardized), but considering
        a 2x1x1 supercell
        """
        import numpy as np         
        
        cell = np.array([(1.,0.,0.),(0.,1.5,0.), (0.1,0.2,2.)])

        positions = ((0.,0.,0.),)
        atomic_numbers = [1]

        # This is the transformation of the cell above
        # (it's already standardized, so it is the identity)
        basic_transf = np.diag([1,1,1])

        self._test_vector_permutations(
            system = (cell, positions, atomic_numbers),
            basic_transf=basic_transf
            )


class TestPaths3D_HPKOT_Supercell(unittest.TestCase):
    """
    Test what happens for a supercell
    """

    def test_supercell(self):
        """
        Test a supercell (BCC).
        This is just a very basic test.
        """
        from seekpath import hpkot

        cell = [[4., 0., 0.], [0., 10., 0.], [0., 0., 4.]]
        positions = [[0., 0., 0.], [0.5, 0.25, 0.5],
                     [0., 0.5, 0.], [0.5, 0.75, 0.5]]
        atomic_numbers = [6, 6, 6, 6]

        system = (cell, positions, atomic_numbers)
        res = hpkot.get_path(system, with_time_reversal=False)

        # Just some basic checks...
        self.assertEqual(res['volume_original_wrt_conv'], 2)
        self.assertEqual(res['volume_original_wrt_prim'], 4)


class TestPaths3D_HPKOT_EdgeCases(unittest.TestCase):
    """
    Test the warnings issued for edge cases
    """

    def basic_test(self, cell, positions, atomic_numbers,
                   check_bravais_lattice, check_string=None):
        """
        Given a cell, the positions and the atomic numbers, checks that
        (only one) warning is issued, of type hpkot.EdgeCaseWarning,
        that the bravais lattice is the expected one,
        and that (optionally, if specified) the warning message contains
        the given string 'check_string'.

        :param cell: 3x3 list of lattice vectors
        :param positions: Nx3 list of (scaled) positions
        :param atomic_number: list of length N with the atomic numbers
        :check_bravais_lattice: a string with the expected Bravais lattice 
            (e.g., 'tI', 'oF', ...)
        :check_string: if specified, this should be contained in the warning
            message
        """
        import warnings

        from seekpath import hpkot

        system = (cell, positions, atomic_numbers)

        with warnings.catch_warnings(record=True) as w:
            res = hpkot.get_path(system, with_time_reversal=False)
            # Checks
            self.assertEqual(res['bravais_lattice'], check_bravais_lattice)
            # Checks on issued warnings
            relevant_w = [_ for _ in w if 
                          issubclass(_.category, hpkot.EdgeCaseWarning)]
            self.assertEqual(len(relevant_w), 1, 
                             'Wrong number of warnings issued! '
                             '({} instead of 1)'.format(len(relevant_w)))
            if check_string is not None:
                self.assertIn(check_string, str(relevant_w[0].message))

    def test_tI(self):
        """
        Test the edge case for tI.
        """
        cell = [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
        positions = [[0., 0., 0.], [0.5, 0.5, 0.5],
                     [0.0, 0.0, 0.1], [0.5, 0.5, 0.6]]
        atomic_numbers = [6, 6, 8, 8]

        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='tI')

    def test_oF_first(self):
        from math import sqrt

        cell = [[sqrt(1. / (1 / 16. + 1 / 25.)), 0., 0.],
                [0., 4., 0.], [0., 0., 5.]]
        positions = [[0., 0., 0.], [0., 0.5, 0.5],
                     [0.5, 0., 0.5], [0.5, 0.5, 0.]]
        atomic_numbers = [6, 6, 6, 6]
        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='oF', check_string="but 1/a^2")

    def test_oF_second(self):
        from math import sqrt

        cell = [[10, 0, 0], [0, 21, 0], [
            0, 0, sqrt(1. / (1 / 100. + 1 / 441.))]]
        positions = [
            [0.1729328200000002,  0.5632488700000001,  0.9531259500000002],
            [0.8270671799999998,  0.4367511299999999,  0.9531259500000002],
            [0.0770671799999998,  0.3132488700000001,  0.7031259500000002],
            [0.9229328200000002,  0.6867511299999999,  0.7031259500000002],
            [0.1729328200000002,  0.0632488700000001,  0.4531259500000002],
            [0.8270671799999998,  0.9367511299999998,  0.4531259500000002],
            [0.0770671799999998,  0.8132488700000001,  0.2031259500000002],
            [0.9229328200000002,  0.1867511299999999,  0.2031259500000002],
            [0.6729328200000002,  0.5632488700000001,  0.4531259500000002],
            [0.3270671799999998,  0.4367511299999999,  0.4531259500000002],
            [0.5770671799999998,  0.3132488700000001,  0.2031259500000002],
            [0.4229328200000002,  0.6867511299999999,  0.2031259500000002],
            [0.6729328200000002,  0.0632488700000001,  0.9531259500000002],
            [0.3270671799999998,  0.9367511299999998,  0.9531259500000002],
            [0.5770671799999998,  0.8132488700000001,  0.7031259500000002],
            [0.4229328200000002,  0.1867511299999999,  0.7031259500000002],
            [0.0000000000000000,  0.5000000000000000,  0.4701481000000003],
            [0.7500000000000000,  0.7500000000000000,  0.2201481000000003],
            [0.0000000000000000,  0.0000000000000000,  0.9701481000000002],
            [0.7500000000000000,  0.2500000000000000,  0.7201481000000003],
            [0.5000000000000000,  0.5000000000000000,  0.9701481000000002],
            [0.2500000000000000,  0.7500000000000000,  0.7201481000000003],
            [0.5000000000000000,  0.0000000000000000,  0.4701481000000003],
            [0.2500000000000000,  0.2500000000000000,  0.2201481000000003]]
        atomic_numbers = [6] * 16 + [8] * 8
        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='oF', check_string="but 1/c^2")

    def test_oI_bc(self):

        cell = [[4., 0., 0.], [0., 5., 0.], [0., 0., 5.]]
        positions = [[0., 0., 0.], [0.5, 0.5, 0.5],
                     [0., 0., 0.1], [0.5, 0.5, 0.6]]
        atomic_numbers = [6, 6, 8, 8]
        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='oI',
                        check_string="but the two longest vectors b and c")

    def test_oC(self):

        cell = [[3., 0., 0.], [0., 3., 0.], [0., 0., 5.]]
        positions = [
            [0.5000000000000000,  0.1136209299999999,  0.7500967299999999],
            [0.5000000000000000,  0.8863790700000000,  0.2500967299999999],
            [0.0000000000000000,  0.6136209300000000,  0.7500967299999999],
            [0.0000000000000000,  0.3863790700000001,  0.2500967299999999],
            [0.0000000000000000,  0.8444605049999999,  0.7659032699999999],
            [0.0000000000000000,  0.1555394950000001,  0.2659032699999999],
            [0.5000000000000000,  0.3444605049999999,  0.7659032699999999],
            [0.5000000000000000,  0.6555394950000001,  0.2659032699999999]]
        atomic_numbers = [6, 6, 6, 6, 8, 8, 8, 8]
        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='oC')

    def test_oA(self):

        cell = [[9., 0., 0.], [0., 3., 0.], [0., 0., 3.]]
        positions = [
            [0.0000000000000000,  0.0000000000000000,  0.0309652399999998],
            [0.0000000000000000,  0.5000000000000000,  0.5309652399999998],
            [0.0000000000000000,  0.5000000000000000,  0.8489601849999999],
            [0.5000000000000000,  0.5000000000000000,  0.1263269549999999],
            [0.0000000000000000,  0.0000000000000000,  0.3489601849999999],
            [0.5000000000000000,  0.0000000000000000,  0.6263269549999999]]
        atomic_numbers = [6, 6, 8, 8, 8, 8]
        self.basic_test(cell, positions, atomic_numbers,
                        check_bravais_lattice='oA')

    # For full coverage, we should also implement the tests for the warnings
    # in the mC lattices, and for the oP warnings.


class TestPaths3D_HPKOT(unittest.TestCase):
    """
    Class to test the creation of paths for all cases using example structures
    """
    # If True, print on stdout the band paths
    verbose_tests = False

    def base_test(self, ext_bravais, with_inv):
        """
        Test a specific extended Bravais symol, 
        with or without inversion (uses the cell whose
        POSCAR is stored in the directories - they have been obtained by 
        Y. Hinuma from the Materials Project).

        :param ext_bravais: a string with the extended Bravais Lattice symbol 
           (like 'cF1', for instance)
        :param with_inv: if True, consider a system with inversion symmetry,
            otherwise one without (in which case, the get_path function is
            called with the kwarg 'with_time_reversal = False')
        """
        import os

        from seekpath import hpkot

        # Get the POSCAR with the example structure
        this_folder = os.path.split(os.path.abspath(hpkot.__file__))[0]
        folder = os.path.join(this_folder, "band_path_data", ext_bravais)
        poscar_with_inv = os.path.join(folder, 'POSCAR_inversion')
        poscar_no_inv = os.path.join(folder, 'POSCAR_noinversion')

        poscar = poscar_with_inv if with_inv else poscar_no_inv
        #asecell = ase.io.read(poscar)
        # system = (asecell.get_cell(), asecell.get_scaled_positions(),
        #    asecell.get_atomic_numbers())
        system = simple_read_poscar(poscar)

        res = hpkot.get_path(system, with_time_reversal=False)

        self.assertEqual(res['bravais_lattice_extended'], ext_bravais)
        self.assertEqual(res['has_inversion_symmetry'], with_inv)

        if self.verbose_tests:
            print("*** {} (inv={})".format(
                ext_bravais, with_inv))
            for p1, p2 in res['path']:
                print("   {} -- {}: {} -- {}".format(p1, p2,
                                                     res['point_coords'][p1], res['point_coords'][p2]))

    def test_aP2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) aP2.
        """
        self.base_test(ext_bravais="aP2", with_inv=True)

    def test_aP2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) aP2.
        """
        self.base_test(ext_bravais="aP2", with_inv=False)

    def test_aP3Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) aP3.
        """
        self.base_test(ext_bravais="aP3", with_inv=True)

    def test_aP3N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) aP3.
        """
        self.base_test(ext_bravais="aP3", with_inv=False)

    def test_cF1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) cF1.
        """
        self.base_test(ext_bravais="cF1", with_inv=True)

    def test_cF1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) cF1.
        """
        self.base_test(ext_bravais="cF1", with_inv=False)

    def test_cF2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) cF2.
        """
        self.base_test(ext_bravais="cF2", with_inv=True)

    def test_cF2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) cF2.
        """
        self.base_test(ext_bravais="cF2", with_inv=False)

    def test_cI1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) cI1.
        """
        self.base_test(ext_bravais="cI1", with_inv=True)

    def test_cI1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) cI1.
        """
        self.base_test(ext_bravais="cI1", with_inv=False)

    def test_cP1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) cP1.
        """
        self.base_test(ext_bravais="cP1", with_inv=True)

    def test_cP1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) cP1.
        """
        self.base_test(ext_bravais="cP1", with_inv=False)

    def test_cP2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) cP2.
        """
        self.base_test(ext_bravais="cP2", with_inv=True)

    def test_cP2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) cP2.
        """
        self.base_test(ext_bravais="cP2", with_inv=False)

    def test_hP1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) hP1.
        """
        self.base_test(ext_bravais="hP1", with_inv=True)

    def test_hP1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) hP1.
        """
        self.base_test(ext_bravais="hP1", with_inv=False)

    def test_hP2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) hP2.
        """
        self.base_test(ext_bravais="hP2", with_inv=True)

    def test_hP2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) hP2.
        """
        self.base_test(ext_bravais="hP2", with_inv=False)

    def test_hR1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) hR1.
        """
        self.base_test(ext_bravais="hR1", with_inv=True)

    def test_hR1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) hR1.
        """
        self.base_test(ext_bravais="hR1", with_inv=False)

    def test_hR2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) hR2.
        """
        self.base_test(ext_bravais="hR2", with_inv=True)

    def test_hR2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) hR2.
        """
        self.base_test(ext_bravais="hR2", with_inv=False)

    def test_mC1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) mC1.
        """
        self.base_test(ext_bravais="mC1", with_inv=True)

    def test_mC1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) mC1.
        """
        self.base_test(ext_bravais="mC1", with_inv=False)

    def test_mC2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) mC2.
        """
        self.base_test(ext_bravais="mC2", with_inv=True)

    def test_mC2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) mC2.
        """
        self.base_test(ext_bravais="mC2", with_inv=False)

    def test_mC3Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) mC3.
        """
        self.base_test(ext_bravais="mC3", with_inv=True)

    def test_mC3N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) mC3.
        """
        self.base_test(ext_bravais="mC3", with_inv=False)

    def test_mP1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) mP1.
        """
        self.base_test(ext_bravais="mP1", with_inv=True)

    def test_mP1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) mP1.
        """
        self.base_test(ext_bravais="mP1", with_inv=False)

# oA1Y does not exist by symmetry
#    def test_oA1Y(self):
#        """
#        Obtain the k-path for a test system with inversion symmetry and
#        Bravais lattice (extended) oA1.
#        """
#        self.base_test(ext_bravais="oA1", with_inv = True)

    def test_oA1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oA1.
        """
        self.base_test(ext_bravais="oA1", with_inv=False)

# oA2Y does not exist by symmetry
#    def test_oA2Y(self):
#        """
#        Obtain the k-path for a test system with inversion symmetry and
#        Bravais lattice (extended) oA2.
#        """
#        self.base_test(ext_bravais="oA2", with_inv = True)

    def test_oA2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oA2.
        """
        self.base_test(ext_bravais="oA2", with_inv=False)

    def test_oC1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oC1.
        """
        self.base_test(ext_bravais="oC1", with_inv=True)

    def test_oC1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oC1.
        """
        self.base_test(ext_bravais="oC1", with_inv=False)

    def test_oC2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oC2.
        """
        self.base_test(ext_bravais="oC2", with_inv=True)

    def test_oC2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oC2.
        """
        self.base_test(ext_bravais="oC2", with_inv=False)

    def test_oF1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oF1.
        """
        self.base_test(ext_bravais="oF1", with_inv=True)

    def test_oF1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oF1.
        """
        self.base_test(ext_bravais="oF1", with_inv=False)

# oF2Y does not exist by symmetry
#    def test_oF2Y(self):
#        """
#        Obtain the k-path for a test system with inversion symmetry and
#        Bravais lattice (extended) oF2.
#        """
#        self.base_test(ext_bravais="oF2", with_inv = True)

    def test_oF2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oF2.
        """
        self.base_test(ext_bravais="oF2", with_inv=False)

    def test_oF3Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oF3.
        """
        self.base_test(ext_bravais="oF3", with_inv=True)

    def test_oF3N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oF3.
        """
        self.base_test(ext_bravais="oF3", with_inv=False)

    def test_oI1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oI1.
        """
        self.base_test(ext_bravais="oI1", with_inv=True)

    def test_oI1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oI1.
        """
        self.base_test(ext_bravais="oI1", with_inv=False)

    def test_oI2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oI2.
        """
        self.base_test(ext_bravais="oI2", with_inv=True)

    def test_oI2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oI2.
        """
        self.base_test(ext_bravais="oI2", with_inv=False)

    def test_oI3Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oI3.
        """
        self.base_test(ext_bravais="oI3", with_inv=True)

    def test_oI3N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oI3.
        """
        self.base_test(ext_bravais="oI3", with_inv=False)

    def test_oP1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) oP1.
        """
        self.base_test(ext_bravais="oP1", with_inv=True)

    def test_oP1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) oP1.
        """
        self.base_test(ext_bravais="oP1", with_inv=False)

    def test_tI1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) tI1.
        """
        self.base_test(ext_bravais="tI1", with_inv=True)

    def test_tI1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) tI1.
        """
        self.base_test(ext_bravais="tI1", with_inv=False)

    def test_tI2Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) tI2.
        """
        self.base_test(ext_bravais="tI2", with_inv=True)

    def test_tI2N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) tI2.
        """
        self.base_test(ext_bravais="tI2", with_inv=False)

    def test_tP1Y(self):
        """
        Obtain the k-path for a test system with inversion symmetry and
        Bravais lattice (extended) tP1.
        """
        self.base_test(ext_bravais="tP1", with_inv=True)

    def test_tP1N(self):
        """
        Obtain the k-path for a test system without inversion symmetry and
        Bravais lattice (extended) tP1.
        """
        self.base_test(ext_bravais="tP1", with_inv=False)

