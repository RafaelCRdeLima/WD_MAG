#!/bin/bash
# Monitor do run norot192 no scicom.
#
# Roda por cron a cada 30 min. Lê o amr_diag.out no destino, mede o dt, projeta
# o término e avisa se algo mudar para pior. Não escreve nada no cluster.
#
# Por que 30 min: o dt caiu ao longo de ~6 h entre t=18 s e t=24 s. Meia hora
# vê a queda com folga e gera 48 linhas de log por dia, não 1440.
#
# Instalação:  crontab -e  ->  */30 * * * * /home/rafael/wd-magnetizada/cluster/scicom/monitor_norot.sh
# Ver o log:   column -s, -t < ~/wd-magnetizada/cluster/scicom/monitor_norot.csv
# Ver alertas: cat ~/wd-magnetizada/cluster/scicom/ALERTAS.txt

set -uo pipefail

HOST=scicom
JOBNAME=rafael-wdmag-norot
RUNDIR='$HOME/wd-mag/runs/dir_norot192'
STOP_TIME=60.0            # stop_time do inputs.norot192
NSTEPS=200                # janela de passos para a média do dt

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CSV="$HERE/monitor_norot.csv"
ALERTAS="$HERE/ALERTAS.txt"
ESTADO="$HERE/.monitor_estado"

# notify-send precisa disto quando chamado pelo cron
export DISPLAY="${DISPLAY:-:0}"
export DBUS_SESSION_BUS_ADDRESS="${DBUS_SESSION_BUS_ADDRESS:-unix:path=/run/user/$(id -u)/bus}"

avisar() {   # avisar <urgencia> <titulo> <corpo>
    local urg="$1" titulo="$2" corpo="$3"
    printf '%s  [%s] %s — %s\n' "$(date '+%F %T')" "$urg" "$titulo" "$corpo" >> "$ALERTAS"
    notify-send -u "$urg" -a "wd-magnetizada" "$titulo" "$corpo" 2>/dev/null || true
}

# ---------------------------------------------------------------- coleta
# Uma única sessão SSH: estado da fila e a cauda do diagnóstico.
DADOS=$(timeout 120 ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" "
    linha=\$(squeue -h -n $JOBNAME -u \$USER -o '%i %T %M %L' 2>/dev/null | head -1)
    echo \"FILA \$linha\"
    if [ -n \"\$linha\" ]; then
        jid=\$(echo \"\$linha\" | awk '{print \$1}')
        echo \"FIM \$(scontrol show job \$jid 2>/dev/null | grep -oP 'EndTime=\K[^ ]+')\"
    fi
    d=$RUNDIR/amr_diag.out
    if [ -f \"\$d\" ]; then
        echo \"MTIME \$(stat -c %Y \"\$d\")\"
        echo \"AGORA \$(date +%s)\"
        echo 'DIAG'
        grep -v '^#' \"\$d\" | tail -$NSTEPS
    fi
" 2>/dev/null)

if [ -z "$DADOS" ]; then
    avisar critical "scicom inacessível" "O monitor não conseguiu ler o run. Rede, fila ou chave."
    printf '%s,ERRO,,,,,,\n' "$(date '+%F %T')" >> "$CSV"
    exit 1
fi

# ---------------------------------------------------------------- análise
LIDO=$(DADOS="$DADOS" python3 - "$STOP_TIME" "$ESTADO" <<'PY'
import sys, datetime

stop_time = float(sys.argv[1])
estado_path = sys.argv[2]
import os
txt = os.environ.get('DADOS', '').split('\n')

fila = fim = ''
mtime = agora = 0
diag = []
modo = None
for l in txt:
    if l.startswith('FILA '): fila = l[5:].strip()
    elif l.startswith('FIM '): fim = l[4:].strip()
    elif l.startswith('MTIME '): mtime = int(l[6:])
    elif l.startswith('AGORA '): agora = int(l[6:])
    elif l.strip() == 'DIAG': modo = 'diag'
    elif modo == 'diag' and l.strip():
        p = l.split()
        if len(p) >= 6:
            try: diag.append((int(p[0]), float(p[1]), float(p[2]), float(p[5])))
            except ValueError: pass

if not diag:
    print('SEMDIAG'); raise SystemExit

step, t, _, _ = diag[-1]
dt_med = sum(d[2] for d in diag) / len(diag)
wt_med = sum(d[3] for d in diag) / len(diag)     # walltime por passo grosso
parado = agora - mtime                           # s desde a última escrita

# taxa e projeção
taxa = (dt_med / wt_med * 3600) if wt_med > 0 else 0   # s_sim por hora
restam = stop_time - t
horas = restam / taxa if taxa > 0 else float('inf')
eta = datetime.datetime.now() + datetime.timedelta(hours=horas) if horas != float('inf') else None

# folga contra o fim da janela de fila
folga = ''
if fim and eta:
    try:
        lim = datetime.datetime.strptime(fim, '%Y-%m-%dT%H:%M:%S')
        folga = f'{(lim - eta).total_seconds()/3600:.0f}'
    except ValueError:
        pass

# dt do check anterior
try:
    dt_ant = float(open(estado_path).read().strip())
except Exception:
    dt_ant = 0.0
open(estado_path, 'w').write(str(dt_med))

print('OK')
print(f'{step}|{t:.4f}|{dt_med:.6e}|{dt_ant:.6e}|{taxa:.4f}|{horas:.1f}|'
      f'{eta:%Y-%m-%d %H:%M}' if eta else f'{step}|{t:.4f}|{dt_med:.6e}|{dt_ant:.6e}|{taxa:.4f}|inf|-')
print(f'{folga}|{parado}|{fila}')
PY
)

if [ "$(head -1 <<<"$LIDO")" = "SEMDIAG" ]; then
    avisar critical "diagnóstico ausente" "amr_diag.out não foi encontrado no run."
    exit 1
fi

IFS='|' read -r STEP T DTMED DTANT TAXA HORAS ETA <<<"$(sed -n 2p <<<"$LIDO")"
IFS='|' read -r FOLGA PARADO FILA                 <<<"$(sed -n 3p <<<"$LIDO")"

ESTADO_FILA=$(awk '{print $2}' <<<"$FILA")

# ---------------------------------------------------------------- registro
[ -f "$CSV" ] || echo "quando,step,t_sim,dt_medio,taxa_s_por_h,horas_restantes,eta,folga_h,estado_fila" > "$CSV"
printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$(date '+%F %T')" "$STEP" "$T" "$DTMED" "$TAXA" "$HORAS" "$ETA" "$FOLGA" "${ESTADO_FILA:-AUSENTE}" >> "$CSV"

# ---------------------------------------------------------------- alertas
# 1. O job saiu da fila. Pode ser sucesso ou morte; quem decide é o t atingido.
if [ -z "$ESTADO_FILA" ]; then
    if python3 -c "import sys; sys.exit(0 if float('$T') >= $STOP_TIME - 0.01 else 1)"; then
        avisar normal "norot192 TERMINOU" "Chegou a t = ${T} s. Alvo era ${STOP_TIME} s."
    else
        avisar critical "norot192 SUMIU DA FILA" "Parou em t = ${T} s de ${STOP_TIME} s. Ver logs/."
    fi
    exit 0
fi

# 2. Escrita parada: mais de 40 min sem tocar o amr_diag (passo grosso ~9 s).
if [ "${PARADO:-0}" -gt 2400 ]; then
    avisar critical "norot192 travado" "amr_diag.out sem escrita há $((PARADO/60)) min, mas o job consta RUNNING."
fi

# 3. Novo colapso de dt: caiu a menos de 60% do check anterior.
if python3 -c "
import sys
a, b = float('$DTANT'), float('$DTMED')
sys.exit(0 if a > 0 and b < 0.60*a else 1)"; then
    avisar critical "dt colapsou de novo" \
        "dt caiu de ${DTANT} para ${DTMED} em t = ${T} s. Taxa agora ${TAXA} s/h, ETA ${ETA}."
fi

# 4. A projeção passou do fim da janela de fila.
if [ -n "$FOLGA" ] && python3 -c "import sys; sys.exit(0 if float('$FOLGA') < 12 else 1)"; then
    avisar critical "norot192 não deve terminar" \
        "ETA ${ETA}, folga de apenas ${FOLGA} h até o fim da janela. t = ${T} s de ${STOP_TIME} s."
fi

exit 0
