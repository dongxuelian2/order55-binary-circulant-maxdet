from fractions import Fraction

from maxdet.mu3 import inner_product_values
from order15_mu3.scripts.verify_color_13_energy9 import minimum_cross_losses


def test_one_edge_thirteen_class_cross_losses() -> None:
    losses = minimum_cross_losses(inner_product_values(15))
    assert losses[0] == Fraction(38, 15)
    assert losses[1] == Fraction(123, 40)
