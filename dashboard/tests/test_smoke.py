"""Teste de fumaca (criterio de aceite do dashboard): com os parametros
padrao e k0=0, a Aba 1 (pages/1_equilibrio.py) deve reproduzir o limite de
Chandrasekhar dentro de 1%. Roda a pagina real via streamlit.testing —
nao e' uma verificacao visual."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

PAGE = Path(__file__).resolve().parent.parent / "pages" / "1_equilibrio.py"


def test_default_params_reproduce_chandrasekhar():
    at = AppTest.from_file(str(PAGE), default_timeout=180)
    at.run()

    assert not at.exception, f"pagina lancou excecao: {at.exception}"

    tables = at.table
    assert len(tables) >= 1, "tabela de escalares nao encontrada"
    scalars = dict(zip(tables[0].value["quantidade"], tables[0].value["valor"]))

    M_msun = float(scalars["M/M☉"])
    rel_err = abs(M_msun - 1.44) / 1.44
    print(f"M = {M_msun:.4f} Msun, erro relativo ao limite de Chandrasekhar: {rel_err:.3%}")
    assert rel_err < 0.01, f"M={M_msun:.4f} Msun, erro {rel_err:.3%} acima de 1%"

    VE = float(scalars["VE"])
    assert VE < 1e-3, f"VE={VE:.3e} acima do V3 do plano"


if __name__ == "__main__":
    test_default_params_reproduce_chandrasekhar()
    print("OK")
