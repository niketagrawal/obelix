import os
import pandas as pd
import pytest
import numpy as np
import shutil
from obelix.descriptor_calculator import Descriptors


@pytest.fixture
def setup_descriptors(request):
    metal_adduct, output_type = request.param
    path_to_input = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data",
        metal_adduct,
        output_type,
    )
    descriptors = Descriptors(
        central_atom="Rh",
        path_to_workflow=path_to_input,
        output_type=output_type,
    )

    if output_type in ["xyz", "crest"]:
        descriptors.calculate_morfeus_descriptors(
            geom_type="BD", solvent=None, printout=False, metal_adduct=metal_adduct
        )
    elif output_type == "gaussian":
        descriptors.calculate_dft_descriptors_from_log(
            geom_type="BD",
            solvent=None,
            extract_xyz_from_log=True,
            printout=False,
            metal_adduct=metal_adduct,
            plot_steric_map=False,
        )
    else:
        raise ValueError("The output_type must be one of 'xyz', 'crest' or 'gaussian'.")

    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "output")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file_path = os.path.join(output_dir, "descriptors.csv")
    descriptors.descriptor_df.to_csv(output_file_path, index=False)

    yield descriptors, metal_adduct, output_type

    shutil.rmtree(output_dir)


@pytest.mark.parametrize(
    "setup_descriptors",
    [
        ("nbd", "xyz"),
        ("pristine", "xyz"),
        ("nbd", "gaussian"),
    ],
    indirect=True,
    ids=[
        "nbd-xyz",
        "pristine-xyz",
        "nbd-gaussian",
    ],  # these ids prints the metal_adduct and output_type in the test name in the pytest report, making it easier to see the pass/fail status of each test case for each combination of metal_adduct and output_type
)
class TestDescriptorCalculation:
    def test_descriptor_values(self, setup_descriptors):
        descriptors, metal_adduct, output_type = setup_descriptors
        self.compare_csv_contents(
            descriptors,
            metal_adduct,
            output_type,
            exclude_columns=["index", "element", "filename_tud"],
        )

    def test_index_values(self, setup_descriptors):
        descriptors, metal_adduct, output_type = setup_descriptors
        self.compare_csv_contents(
            descriptors, metal_adduct, output_type, include_columns=["index"]
        )

    def test_element_values(self, setup_descriptors):
        descriptors, metal_adduct, output_type = setup_descriptors
        self.compare_csv_contents(
            descriptors, metal_adduct, output_type, include_columns=["element"]
        )

    def test_filename_values(self, setup_descriptors):
        descriptors, metal_adduct, output_type = setup_descriptors
        self.compare_csv_contents(
            descriptors, metal_adduct, output_type, include_columns=["filename_tud"]
        )

    def compare_csv_contents(
        self,
        descriptors,
        metal_adduct,
        output_type,
        include_columns=None,
        exclude_columns=None,
    ):
        """
        This is a helper function to compare the contents of the output csv file with the expected csv file. Extracting out this function helps to avoid code duplication in the test methods.
        """
        path_expected_csv = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "expected_output",
            f"descriptors_{metal_adduct}_{output_type}.csv",
        )
        output_csv = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "output",
            "descriptors.csv",
        )

        expected_df = pd.read_csv(path_expected_csv)
        output_df = pd.read_csv(output_csv)

        # filter the columns based on the include_columns or exclude_columns parameter and compare the values. If include_columns is provided, use DataFrame.equals for comparison. If exclude_columns is provided, use np.allclose for comparison. This is because DataFrame.equals is used for non-numeric data comparison (index, filenme and element values), while np.allclose is used for numeric data comparison (descriptor values).
        if include_columns:
            filtered_columns = expected_df.columns[
                expected_df.columns.str.contains("|".join(include_columns))
            ]
            expected_df = expected_df.loc[:, filtered_columns]
            output_df = output_df.loc[:, filtered_columns]
            # Use DataFrame.equals for non-numeric data comparison
            assert expected_df.equals(
                output_df
            ), f"The {include_columns} values in the output csv file do not match the expected values for {metal_adduct} with {output_type}."
        elif exclude_columns:
            # For numeric data, use np.allclose after excluding non-numeric columns
            filtered_columns = expected_df.columns[
                ~expected_df.columns.str.contains("|".join(exclude_columns))
            ]
            expected_df = expected_df.loc[:, filtered_columns]
            output_df = output_df.loc[:, filtered_columns]

            output_descriptor_values = output_df.to_numpy()

            # Convert True/False values to 1/0
            output_descriptor_values = np.where(
                output_descriptor_values == True, 1, output_descriptor_values
            )
            output_descriptor_values = np.where(
                output_descriptor_values == False, 0, output_descriptor_values
            )

            # Convert each element in the array to a numeric type (float by default). If an element cannot be converted (e.g., a string that doesn't represent a number), it is replaced with NaN (Not a Number), as a result of the errors='coerce' parameter. This step ensures that the entire array is numeric, which is a prerequisite for using np.allclose.
            output_descriptor_values = pd.to_numeric(
                output_descriptor_values.flatten(), errors="coerce"
            ).reshape(output_descriptor_values.shape)

            expected_descriptor_values = expected_df.to_numpy()

            expected_descriptor_values = np.where(
                expected_descriptor_values == True, 1, expected_descriptor_values
            )
            expected_descriptor_values = np.where(
                expected_descriptor_values == False, 0, expected_descriptor_values
            )

            expected_descriptor_values = pd.to_numeric(
                expected_descriptor_values.flatten(), errors="coerce"
            ).reshape(expected_descriptor_values.shape)

            assert np.allclose(
                output_descriptor_values, expected_descriptor_values, equal_nan=True
            ), f"The descriptor values in the output csv file do not match the expected descriptor values for {metal_adduct} with {output_type}."
        else:
            raise ValueError(
                "Either include_columns or exclude_columns must be provided."
            )
