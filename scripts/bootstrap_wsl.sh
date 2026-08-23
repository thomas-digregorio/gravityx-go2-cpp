#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "${repo_root,,}" in
  *onedrive*) echo "refusing OneDrive repository path: ${repo_root}" >&2; exit 2 ;;
esac

tool_root="${HOME}/.local/share/gravityx-go2-cpp"
env_root="${tool_root}/env"
gravity_root="${repo_root}/.deps/Gravity"
gravity_revision="f5af33ed9572829d96c6f54b1a3ad30f53677fe7"

mkdir -p "${tool_root}/bin" "${tool_root}/cache" "${repo_root}/.deps"
if [[ ! -x "${tool_root}/bin/micromamba" ]]; then
  curl -L --fail --retry 3 https://micro.mamba.pm/api/micromamba/linux-64/latest \
    -o "${tool_root}/cache/micromamba.tar.bz2"
  tar -xjf "${tool_root}/cache/micromamba.tar.bz2" -C "${tool_root}" bin/micromamba
fi
if [[ ! -x "${env_root}/bin/cmake" ]]; then
  MAMBA_ROOT_PREFIX="${tool_root}/mamba-root" "${tool_root}/bin/micromamba" create -y \
    -p "${env_root}" -c conda-forge cmake ninja ipopt pkg-config nlohmann_json
fi
if [[ ! -d "${gravity_root}/.git" ]]; then
  git clone https://github.com/coin-or/Gravity.git "${gravity_root}"
fi
git -C "${gravity_root}" fetch origin "${gravity_revision}"
git -C "${gravity_root}" checkout --detach "${gravity_revision}"

export PATH="${env_root}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export IPOPT_ROOT_DIR="${env_root}"
export LD_LIBRARY_PATH="${env_root}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

cmake --fresh -S "${gravity_root}" -B "${gravity_root}/build-codex" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DIpopt=ON -DMP=OFF -DOPT_PARSER=OFF -DGurobi=OFF -DCplex=OFF \
  -DBonmin=OFF -DClp=OFF -DMosek=OFF -DSdpa=OFF -DXlnt=OFF \
  -DBoost=OFF -DQpp=OFF -DOpenMPI=OFF \
  -DCMAKE_PREFIX_PATH="${env_root}"
cmake --build "${gravity_root}/build-codex" --target gravity -j "$(nproc)"

cmake --fresh -S "${repo_root}" -B "${repo_root}/.build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${env_root}" \
  -DIPOPT_ROOT="${env_root}"
cmake --build "${repo_root}/.build" -j "$(nproc)"
ctest --test-dir "${repo_root}/.build" --output-on-failure
