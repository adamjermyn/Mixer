~/.local/bin/gprof2dot -s Build/out > profile.dot
dot -Tps profile.dot -o profile.eps

