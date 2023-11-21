import pygame
from dataclasses import dataclass
import numpy as np
from utils import ScenarioType
import do_mpc


@dataclass
class BioGasPlantVis:
    max_volume: float
    _screen: pygame.surface.Surface

    _color = "black"

    def __post_init__(self):
        self._screen_size = self._screen.get_size()

        self._bb_top_left = (self._screen_size[0] * 0.1, self._screen_size[1] * 0.1)
        self._bb_bottom_right = (self._screen_size[0] * 0.5, self._screen_size[1] * 0.9)
        self._bb_bottom_left = (self._bb_top_left[0], self._bb_bottom_right[1])
        self._bb_top_right = (self._bb_bottom_right[0], self._bb_top_left[1])

        self._draw_width = round(self._screen_size[0] / 400)

    def draw(self, v_max, v_ch4, v_co2, v_h2o, u_actual):
        height_gas_storage = (self._bb_bottom_left[1] - self._bb_top_left[1]) / 4.0

        top_left_gas = list(self._bb_top_left)
        top_left_gas[1] += (self._bb_bottom_left[1] - self._bb_top_left[1]) / 3.0

        top_left_liquid = list(self._bb_top_left)
        top_left_liquid[1] += (self._bb_bottom_left[1] - self._bb_top_left[1]) / 2.0

        top_left_feed = list(self._bb_top_left)
        top_left_feed[1] += (self._bb_bottom_left[1] - self._bb_top_left[1]) * 3.0 / 4.0

        width = self._bb_bottom_right[0] - self._bb_bottom_left[0]
        height_gas = top_left_liquid[1] - top_left_gas[1]
        height_liquid = top_left_feed[1] - top_left_liquid[1]
        height_feed = self._bb_bottom_left[1] - top_left_feed[1]

        self._draw_fermenter(
            top_left_gas, width, height_gas, top_left_liquid, height_liquid
        )
        self._draw_gas_storage(width, height_gas_storage, v_max, v_ch4, v_co2, v_h2o)
        self._draw_feed(top_left_feed, width, height_feed, u_actual)

    def _draw_fermenter(
        self, top_left_gas, width, height_gas, top_left_liquid, height_liquid
    ):
        pygame.draw.rect(
            self._screen,
            self._color,
            (*top_left_gas, width, height_gas),
            width=self._draw_width,
        )

        pygame.draw.rect(
            self._screen,
            self._color,
            (*top_left_liquid, width, height_liquid),
            width=self._draw_width,
        )

    def _draw_gas_storage(self, width, height_gas_storage, v_max, v_ch4, v_co2, v_h2o):
        v_total_fill = v_ch4 + v_co2 + v_h2o

        pygame.draw.rect(
            self._screen,
            self._color,
            (*self._bb_top_left, width, height_gas_storage),
            width=self._draw_width,
        )

        # Draw fillings

        top_left = list(self._bb_top_left)
        height_filling = v_total_fill / v_max * height_gas_storage
        top_left[1] += height_gas_storage - height_filling

        # Draw filling with CH4
        width_ch4 = v_ch4 / v_total_fill * width
        pygame.draw.rect(
            self._screen,
            "green",
            (*top_left, width_ch4, height_filling),
            width=0,
        )

        # Draw filling with CO2
        width_co2 = v_co2 / v_total_fill * width
        top_left[0] += width_ch4
        pygame.draw.rect(
            self._screen,
            "black",
            (*top_left, width_co2, height_filling),
            width=0,
        )

        # Draw filling with H2O
        width_h2o = v_h2o / v_total_fill * width
        top_left[0] += width_co2
        pygame.draw.rect(
            self._screen,
            "blue",
            (*top_left, width_h2o, height_filling),
            width=0,
        )

    def _draw_feed(self, top_left_feed, width, height_feed, u):
        width = width / (len(u) * 2 - 1)
        for u in u:
            height = u * height_feed
            top_left_feed[1] += height_feed - height
            pygame.draw.rect(
                self._screen,
                "brown",
                (*top_left_feed, width, height),
                width=0,
            )
            top_left_feed[0] += 2 * width
            top_left_feed[1] -= height_feed - height


@dataclass
class DataVis:
    _screen: pygame.surface.Surface

    _color = "black"

    def __post_init__(self):
        self._screen_size = self._screen.get_size()

        self._top_left = (self._screen_size[0] * 0.6, self._screen_size[1] * 0.1)
        self._bottom_right = (self._screen_size[0] * 0.9, self._screen_size[1] * 0.9)
        self._bottom_left = (self._top_left[0], self._bottom_right[1])
        self._top_right = (self._bottom_right[0], self._top_left[1])

        self._width = self._top_right[0] - self._top_left[0]
        self._height = self._bottom_left[1] - self._top_left[1]

        self._draw_width = round(self._screen_size[0] / 400)

    def draw(self, x: np.ndarray, state_names: list, y: np.ndarray, meas_names: list):
        pygame.draw.rect(
            self._screen,
            self._color,
            (*self._top_left, self._width, self._height),
            width=self._draw_width,
            border_radius=self._draw_width * 3,
        )

        # Create a font object
        font_title = pygame.font.Font(None, round(self._screen_size[1] / 20))
        font_data = pygame.font.Font(None, round(self._screen_size[1] / 30))

        # Render text
        text_title = font_title.render("States", True, "black")
        x = x.tolist()
        text_data = [
            font_data.render(f"{state_name} = {_x[0]:.2f}", True, "black")
            for state_name, _x in zip(state_names, x)
        ]

        # Draw text to screen
        self._screen.blit(
            text_title,
            (self._top_left[0] * 1.02, self._top_left[1] * 1.05),
        )
        for idx, text in enumerate(text_data):
            self._screen.blit(
                text,
                (
                    self._top_left[0] * 1.02,
                    self._top_left[1] * 1.05
                    + (idx + 2) * round(self._screen_size[1] / 30),
                ),
            )

        text_title = font_title.render("Measurements", True, "black")
        # Draw text to screen
        self._screen.blit(
            text_title,
            (
                self._top_left[0] + self._width / 2 * 1.05,
                self._top_left[1] * 1.05,
            ),
        )
        y = y.tolist()
        text_data = [
            font_data.render(f"{meas_names} = {_y:.2f}", True, "black")
            for meas_names, _y in zip(meas_names, y)
        ]

        for idx, text in enumerate(text_data):
            self._screen.blit(
                text,
                (
                    self._top_left[0] + self._width / 2 * 1.05,
                    self._top_left[1] * 1.05
                    + (idx + 2) * round(self._screen_size[1] / 30),
                ),
            )


def visualize(
    bga: BioGasPlantVis,
    data: DataVis,
    state_names: list[str],
    meas_names: list[str],
    Tx: np.ndarray,
    Ty: np.ndarray,
    x0_norm: np.ndarray,
    simulator: do_mpc.simulator.Simulator,
    u_actual: np.ndarray,
    V_GAS_STORAGE_MAX: float,
    scenario_type: ScenarioType,
):
    x0 = np.copy(x0_norm)
    x0 *= np.array([Tx]).T
    y = np.array([simulator.data._aux[-1, idx + 1] * Ty for idx, Ty in enumerate(Ty)])
    if scenario_type == ScenarioType.COGENERATION:
        bga.draw(
            V_GAS_STORAGE_MAX,
            x0[-2][0],
            x0[-1][0],
            simulator.data._aux[-1, 9],
            u_actual.flatten(),
        )
    data.draw(x0, state_names, y, meas_names)

    # flip() the display to put your work on screen
    pygame.display.flip()
