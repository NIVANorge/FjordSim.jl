#!/usr/bin/env bash
#
# Make the VM able to run the FjordSim container. Rendered by `fjordsim-gcp bootstrap` (envsubst),
# scp'd to the VM and run there.
#
# The launcher job creates a plain `ubuntu-2204-lts` instance, so none of this is preinstalled. Ask
# the admin to switch the job to a Deep Learning VM image (see "Asking the admin" in gcp/README.md)
# and every stage below turns into a no-op — which is also what happens on the second call, and
# after `down`/`up`, since the boot disk survives a guest poweroff.
#
# It does not survive a *preemption*: that deletes the instance and its disk, and bootstrap has to
# run again from scratch.
set -euo pipefail

readonly IMAGE_URI="${IMAGE_URI}"
readonly STAGE=/mnt/stage
readonly REGISTRY_HOST="${IMAGE_URI%%/*}"
readonly REBOOT_STAMP="$HOME/.fjordsim-bootstrap-rebooted"

apt_get() { sudo DEBIAN_FRONTEND=noninteractive apt-get -qq "$@"; }

# --- 1. NVIDIA driver ---------------------------------------------------------------------------
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    echo "[1/5] driver: already present"
elif [[ -f "$REBOOT_STAMP" ]]; then
    # We installed it and rebooted once already and the GPU is still not there. Another reboot
    # would just loop.
    echo "[1/5] driver: installed but nvidia-smi still fails after a reboot" >&2
    nvidia-smi || true
    exit 1
else
    echo "[1/5] driver: installing (this is the slow part, ~5 min)"
    apt_get update
    # The driver is a DKMS module, so the toolchain has to be in place *before* the driver package
    # is configured, and it has to be the compiler that built this kernel. Both halves matter: dpkg
    # configures nvidia-dkms-open before the `gcc` metapackage that provides /usr/bin/cc, so a
    # single transaction fails DKMS's CC sanity check; and jammy's default gcc-11 does not know the
    # -ftrivial-auto-var-init=zero that the gcc-12-built -gcp kernel's config demands, so the build
    # fails again even once `cc` exists. Priority 100 beats the gcc-11 entry's 20.
    kernel_gcc="gcc-$(sed -n 's/.*x86_64-linux-gnu-gcc-\([0-9]*\).*/\1/p' /proc/version)"
    [[ "$kernel_gcc" != "gcc-" ]] || { echo "cannot read the kernel's compiler from /proc/version" >&2; exit 1; }
    apt_get install -y build-essential dkms "linux-headers-$(uname -r)" "$kernel_gcc"
    sudo update-alternatives --install /usr/bin/cc cc "/usr/bin/$kernel_gcc" 100
    curl -fsSLO https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
    sudo dpkg -i cuda-keyring_1.1-1_all.deb
    rm -f cuda-keyring_1.1-1_all.deb
    apt_get update
    apt_get install -y cuda-drivers
    touch "$REBOOT_STAMP"
    # Do not reboot from here: this script is running inside the SSH session that the reboot would
    # kill, so its exit status would be lost. 75 is EX_TEMPFAIL — `fjordsim-gcp bootstrap` reads it
    # as "reboot me and call again", which keeps the signal unambiguous.
    echo '[1/5] driver: installed; needs a reboot to load the kernel module'
    exit 75
fi

# --- 2. Docker ----------------------------------------------------------------------------------
if command -v docker >/dev/null 2>&1; then
    echo "[2/5] docker: already present"
else
    echo "[2/5] docker: installing"
    curl -fsSL https://get.docker.com | sudo sh
    sudo usermod -aG docker "$USER"
fi

# --- 3. nvidia-container-toolkit ----------------------------------------------------------------
if command -v nvidia-ctk >/dev/null 2>&1; then
    echo "[3/5] container toolkit: already present"
else
    echo "[3/5] container toolkit: installing"
    curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey \
        | sudo gpg --batch --yes --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
    curl -fsSL https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list \
        | sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
        | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list >/dev/null
    apt_get update
    apt_get install -y nvidia-container-toolkit
    sudo nvidia-ctk runtime configure --runtime=docker
    sudo systemctl restart docker
fi

# --- 4. stage directory and registry auth -------------------------------------------------------
echo "[4/5] stage dir and registry auth"
sudo mkdir -p "$STAGE"/{data,results,configs}
sudo chown -R "$USER" "$STAGE"
gcloud auth configure-docker "$REGISTRY_HOST" --quiet

# --- 5. verify ----------------------------------------------------------------------------------
# `sg docker` picks up the group membership without a logout, which the freshly-added user does not
# have in this session yet.
echo "[5/5] verifying"
nvidia-smi
sg docker -c 'docker info --format "{{.Runtimes}}"' | grep -q nvidia \
    || { echo "docker has no nvidia runtime configured" >&2; exit 1; }

if ! sg docker -c "docker pull '$IMAGE_URI'"; then
    # Whoever we are is a question only the metadata server can answer authoritatively — the
    # launcher job picks the identity, and nothing on this side is told which one it chose.
    identity="$(curl -sf -H 'Metadata-Flavor: Google' \
        http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/email \
        || echo '(could not reach the metadata server)')"
    cat >&2 <<MSG

docker pull failed. If that was a 403, the identity this VM runs as
  $identity
needs roles/artifactregistry.reader on the repository behind
  $IMAGE_URI
which is an admin request — see "Asking the admin" in gcp/README.md. Until then you can move the
image through the bucket instead: \`docker save | gzip\` locally, upload, and \`docker load\` here.
MSG
    exit 1
fi

rm -f "$REBOOT_STAMP"
echo "bootstrap complete"
