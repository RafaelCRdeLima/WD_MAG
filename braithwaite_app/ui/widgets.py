"""Small reusable Qt widgets shared across the app's screens."""

from PyQt6.QtWidgets import QDoubleSpinBox


class ScientificDoubleSpinBox(QDoubleSpinBox):
    """A QDoubleSpinBox that displays and accepts scientific notation
    (e.g. "1.000e+09") instead of the plain fixed-point form Qt uses by
    default -- reading a bare density like 1000000000.000 g/cm3 digit by
    digit is exactly the kind of thing that hides a typo (an extra or
    missing zero) until the wrong star gets built.

    Step behavior is multiplicative (x10 / x0.1) rather than additive --
    additive steps are meaningless across a range spanning many decades
    (1e5 to 1e13 g/cm3).
    """

    def __init__(self, decimals: int = 3, parent=None):
        super().__init__(parent)
        self._decimals = decimals
        self.setDecimals(10)  # internal float precision; display is overridden below

    def textFromValue(self, value: float) -> str:
        if value == 0:
            return "0"
        return f"{value:.{self._decimals}e}"

    def valueFromText(self, text: str) -> float:
        try:
            return float(text.strip().replace(",", "."))
        except ValueError:
            return self.value()

    def stepBy(self, steps: int):
        factor = 10 ** steps
        new_value = self.value() * factor
        new_value = max(self.minimum(), min(self.maximum(), new_value))
        self.setValue(new_value)

    def validate(self, text: str, pos: int):
        # Accept anything that could become a valid float while typing
        # (Qt's default validator rejects "1e" mid-entry, "1.2e-" etc.)
        from PyQt6.QtGui import QValidator
        t = text.strip()
        if t in ("", "-", "+"):
            return (QValidator.State.Intermediate, text, pos)
        try:
            float(t.replace(",", "."))
            return (QValidator.State.Acceptable, text, pos)
        except ValueError:
            partial_chars = set("0123456789.eE+-,")
            if all(c in partial_chars for c in t):
                return (QValidator.State.Intermediate, text, pos)
            return (QValidator.State.Invalid, text, pos)
