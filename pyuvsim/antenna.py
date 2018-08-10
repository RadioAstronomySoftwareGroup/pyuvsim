
class Antenna(object):
    @profile
    def __init__(self, name, number, enu_position, beam_id):
        self.name = name
        self.number = number
        # ENU position in meters relative to the telescope_location
        self.pos_enu = enu_position * units.m
        # index of beam for this antenna from array.beam_list
        self.beam_id = beam_id

    @profile
    def get_beam_jones(self, array, source_az_za, frequency):
        # get_direction_jones needs to be defined on UVBeam
        # 2x2 array of Efield vectors in Az/ZA
        # return array.beam_list[self.beam_id].get_direction_jones(source_lmn, frequency)

        source_az = np.array([source_az_za[0]])
        source_za = np.array([source_az_za[1]])
        freq = np.array([frequency.to('Hz').value])

        if array.beam_list[self.beam_id].data_normalization != 'peak':
            array.beam_list[self.beam_id].peak_normalize()
        array.beam_list[self.beam_id].interpolation_function = 'az_za_simple'

        interp_data, interp_basis_vector = \
            array.beam_list[self.beam_id].interp(az_array=source_az,
                                                 za_array=source_za,
                                                 freq_array=freq)

        # interp_data has shape: (Naxes_vec, Nspws, Nfeeds, 1 (freq), 1 (source position))
        jones_matrix = np.zeros((2, 2), dtype=np.complex)
        # first axis is feed, second axis is theta, phi (opposite order of beam!)
        jones_matrix[0, 0] = interp_data[1, 0, 0, 0, 0]
        jones_matrix[1, 1] = interp_data[0, 0, 1, 0, 0]
        jones_matrix[0, 1] = interp_data[0, 0, 0, 0, 0]
        jones_matrix[1, 0] = interp_data[1, 0, 1, 0, 0]

        return jones_matrix

    def __eq__(self, other):
        return ((self.name == other.name)
                and np.allclose(self.pos_enu.to('m').value, other.pos_enu.to('m').value, atol=1e-3)
                and (self.beam_id == other.beam_id))

    def __gt__(self, other):
        return (self.number > other.number)

    def __ge__(self, other):
        return (self.number >= other.number)

    def __lt__(self, other):
        return not self.__ge__(other)

    def __le__(self, other):
        return not self.__gt__(other)
