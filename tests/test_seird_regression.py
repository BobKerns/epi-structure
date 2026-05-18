import math
import unittest

from epi_structure import (
    DiseaseParameters,
    EpidemicModel,
    PopulationParameters,
    SimulationParameters,
    StructuredEpidemicModel,
)


class TestSEIRDRegression(unittest.TestCase):
    def setUp(self) -> None:
        self.sim = SimulationParameters(duration=40.0, time_step=0.1, output_stride=20)

    def test_single_population_seir_has_no_deceased_column(self) -> None:
        disease = DiseaseParameters(latent_period=5.0)
        pop = PopulationParameters(
            name="A",
            size=10_000,
            beta=0.4,
            initial_infected=10,
            disease=disease,
        )

        df = EpidemicModel(population=pop, simulation=self.sim).simulate(tidy=True)

        self.assertNotIn("deceased", df.columns)
        self.assertIn("total_population", df.columns)
        self.assertGreater(float(df["infected"].max()), 0.0)

    def test_single_population_seird_tracks_deceased_and_conserves_total(self) -> None:
        disease = DiseaseParameters(
            latent_period=5.0,
            compartments=["S", "E", "I", "R", "D"],
            case_fatality_rate=0.01,
        )
        pop = PopulationParameters(
            name="A",
            size=10_000,
            beta=0.4,
            initial_infected=10,
            disease=disease,
        )

        df = EpidemicModel(population=pop, simulation=self.sim).simulate(tidy=True)

        self.assertIn("deceased", df.columns)
        self.assertGreater(float(df["deceased"].iloc[-1]), 0.0)

        total = df["susceptible"] + df["exposed"] + df["infected"] + df["recovered"] + df["deceased"]
        for value in total:
            self.assertTrue(math.isclose(float(value), 10_000.0, rel_tol=0.0, abs_tol=1e-5))

    def test_structured_seird_tracks_deceased(self) -> None:
        disease = DiseaseParameters(
            latent_period=5.0,
            compartments=["S", "E", "I", "R", "D"],
            case_fatality_rate=0.01,
        )
        p1 = PopulationParameters(name="A", size=8_000, beta=0.35, initial_infected=8, disease=disease)
        p2 = PopulationParameters(name="B", size=6_000, beta=0.30, initial_infected=5, disease=disease)

        model = StructuredEpidemicModel(
            populations=[p1, p2],
            contact_matrix=[[0.5, 0.1], [0.1, 0.4]],
            simulation=SimulationParameters(duration=100.0, time_step=0.1, output_stride=50),
        )

        df = model.simulate(tidy=True)

        self.assertIn("deceased", df.columns)

        last = df.sort_values("time").groupby("population").tail(1)
        self.assertGreater(float(last["deceased"].sum()), 0.0)


if __name__ == "__main__":
    unittest.main()
