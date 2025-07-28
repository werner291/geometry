{
  description = "Motion Planning Research";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    pre-commit-hooks = { url = "github:cachix/pre-commit-hooks.nix"; };
    rust-overlay = { url = "github:oxalica/rust-overlay"; };
  };
  outputs = { self, nixpkgs, pre-commit-hooks, rust-overlay }:
    let
      supportedSystems = [ "x86_64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
      pkgsFor = system: import nixpkgs {
        inherit system;
        overlays = [ rust-overlay.overlays.default ];
      };
    in
    {
      formatter.x86_64-linux = nixpkgs.legacyPackages.x86_64-linux.nixpkgs-fmt;
      devShells = forAllSystems (system: {
        default = (pkgsFor system).callPackage ./shell.nix { inherit pre-commit-hooks; inherit system; };
      });
      packages = forAllSystems (system: {
        default = (pkgsFor system).callPackage ./. { };
      });
    };
}
