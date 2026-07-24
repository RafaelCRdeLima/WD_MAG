"""Protocolo comum dos termos de Bernoulli plugaveis em scf.hachisu_scf().

Cada ingrediente de fisica (rotacao, poloidal, toroidal autoconsistente)
contribui um termo aditivo ao MESMO Bernoulli:

    H + Phi - C_rot(varpi) - M_pol(u) + M_tor(rho*varpi^2) = C

rearranjado para o loop (a forma que scf.py usa):

    H = C - Phi + C_rot(varpi) + M_pol(u) - M_tor(rho*varpi^2)

Cada termo abaixo devolve sua contribuicao JA' COM O SINAL CERTO para ser
SOMADA no lado direito acima -- scf.hachisu_scf() so' soma o que nao for
None. None = termo desligado (contribuicao e energia identicamente zero).

Dois tipos de termo:

- EXPLICITOS (rotation, poloidal): o potencial e' um campo na malha que
  NAO depende do rho_novo sendo resolvido na iteracao atual -- rotacao
  depende so' da geometria fixa da malha (varpi); poloidal depende de u,
  resolvido a partir do rho da iteracao ANTERIOR via Grad-Shafranov (o
  mesmo defasamento que ja existia antes desta arquitetura). Interface:

      term.update(rho, r, theta, **kwargs)   # recalcula estado interno
      term.potential(r, theta) -> campo      # contribuicao ao Bernoulli
      term.energy(rho, r, theta) -> dict     # para virial/diagnosticos

- IMPLICITO (toroidal_sc): seu potencial depende do proprio rho sendo
  resolvido, NO MESMO ponto de malha, na MESMA iteracao -- nao ha' defasagem,
  porque o ramo toroidal autoconsistente existe exatamente para o campo
  responder algebricamente a' densidade local. Isso transforma o passo 9 do
  loop de inversao direta da EOS numa busca de raiz por ponto (ver
  scf.py :: _solve_rho_implicit). Interface:

      term.potential_of_rho(rho_trial, r, theta) -> campo  # depende de rho
      term.energy(rho, r, theta) -> dict (uma vez que rho e' conhecido)

poloidal e toroidal sao MUTUAMENTE EXCLUSIVOS -- campo misto
poloidal+toroidal autoconsistente esta' fora de escopo (a barotropia nao
entrega o Bt/Bp desejado; e' por isso que existe D6 em
plano_wd_magnetizada.md, que fica INTACTO: toroidal.py::impose_toroidal()
continua sendo a imposicao a posteriori para o dado inicial do Castro, um
objeto completamente diferente deste pacote). scf.hachisu_scf() garante a
exclusao mutua.
"""
