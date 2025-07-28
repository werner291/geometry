{ pkgs, pre-commit-hooks, system }:
let
  toolchain = pkgs.rust-bin.stable.latest.default.override {
    extensions = [ "rust-src" "rustfmt" "clippy" ];
  };
  pre-commit-check = pre-commit-hooks.lib.${system}.run {
    src = ./.;
    hooks = {
      nixpkgs-fmt.enable = true;
      rustfmt = {
        enable = true;
        package = toolchain;
      };
    };
  };
in
pkgs.mkShell rec {
  buildInputs = with pkgs; [
    toolchain
    xorg.libXcursor
    xorg.libXrandr
    xorg.libXi
    xorg.libX11
    libxkbcommon
    libGL
    alsa-lib
    pkg-config
    cmake
    udev
    vulkan-loader
    rustPlatform.bindgenHook
    ffmpeg
    fontconfig
    (pkgs.python312.withPackages (python-pkgs: with python-pkgs; [
      pandas
      numpy
      matplotlib
      seaborn
      msgpack #
    ]))
    bashInteractive
  ] ++ pre-commit-check.enabledPackages;

  shellHook = pre-commit-check.shellHook;

  LD_LIBRARY_PATH = "${pkgs.lib.makeLibraryPath buildInputs}";
  RUST_SRC_PATH = "${toolchain}/lib/rustlib/src/rust/library";
}
