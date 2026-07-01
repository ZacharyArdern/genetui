#!/usr/bin/env sh
set -e

REPO="ZacharyArdern/genetui"
BIN="genetui"
INSTALL_DIR="${INSTALL_DIR:-$HOME/.local/bin}"

# Detect OS and arch
OS="$(uname -s)"
ARCH="$(uname -m)"

case "${OS}" in
  Darwin*)
    ASSET="genetui-macos-arm64.tar.gz"
    ;;
  Linux*)
    case "${ARCH}" in
      x86_64) ASSET="genetui-linux-x86_64.tar.gz" ;;
      *)
        echo "Unsupported Linux architecture: ${ARCH}" >&2
        echo "Please build from source: https://github.com/${REPO}" >&2
        exit 1
        ;;
    esac
    ;;
  *)
    echo "Unsupported OS: ${OS}" >&2
    echo "Please build from source: https://github.com/${REPO}" >&2
    exit 1
    ;;
esac

# Get latest release tag
TAG="$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" | grep '"tag_name"' | sed 's/.*"tag_name": *"\([^"]*\)".*/\1/')"
if [ -z "${TAG}" ]; then
  echo "Could not determine latest release tag." >&2
  exit 1
fi

URL="https://github.com/${REPO}/releases/download/${TAG}/${ASSET}"

echo "Installing ${BIN} ${TAG} for ${OS}/${ARCH}..."
echo "Downloading ${URL}"

TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

curl -fsSL "${URL}" -o "${TMP}/${ASSET}"
tar -xzf "${TMP}/${ASSET}" -C "${TMP}"

mkdir -p "${INSTALL_DIR}"
mv "${TMP}/${BIN}" "${INSTALL_DIR}/${BIN}"
chmod +x "${INSTALL_DIR}/${BIN}"

echo ""
echo "Installed ${BIN} to ${INSTALL_DIR}/${BIN}"

# Warn if INSTALL_DIR not in PATH
case ":${PATH}:" in
  *":${INSTALL_DIR}:"*) ;;
  *)
    echo ""
    echo "Note: ${INSTALL_DIR} is not in your PATH."
    echo "Add this to your shell profile:"
    echo "  export PATH=\"${INSTALL_DIR}:\$PATH\""
    ;;
esac
