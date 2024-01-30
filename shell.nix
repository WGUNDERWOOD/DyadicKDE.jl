let
pkgs = import <nixpkgs> { };
in pkgs.mkShell {
    buildInputs = with pkgs; [
        julia
        python3
    ];
}