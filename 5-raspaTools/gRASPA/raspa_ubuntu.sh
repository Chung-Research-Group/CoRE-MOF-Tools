export RASPA_DIR=/home/carbon/miniconda3/envs/RASPA2

MAX_JOBS=120

function limit_jobs {
  while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
    sleep 600
  done
}

function run_simulation {
  local struc="$1"
  local P="$2"
  ${RASPA_DIR}/bin/simulate simulation.input
  echo "finished" > "${struc}-${P}Pa"
}

DIR=$(pwd)
mkdir -p calculation
cd calculation

while read struc
do
  mkdir -p "${struc}"
  cd "${struc}"

  UNIT_CELL=$(grep -w "$struc" "${DIR}/mof_list_unitcell.dat" | awk '{print "\t"$2"\t"$3"\t"$4}')


    while read P
    do
      mkdir -p "${P}Pa"
      cd "${P}Pa"
      cp "${DIR}/CoREMOF2019_public_v2_20241118/CR/${struc}.cif" ./
      cp "${DIR}/simulation.input" ./

      sed -i "s/^.*FrameworkName.*$/FrameworkName $struc/" simulation.input
      sed -i "s/^.*UnitCells.*$/UnitCells $UNIT_CELL/" simulation.input
      sed -i "s/^.*ExternalPressure.*$/ExternalPressure $P/" simulation.input

      limit_jobs

      run_simulation "$struc" "$P" &

      cd ..
    done < "${DIR}/p_list"
  cd ..
done < ${DIR}/$1

wait