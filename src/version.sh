#!/bin/bash

GIT_VER=$(git rev-parse --short HEAD)

update_ver() {
	cat src/version.h.in | sed s,{{GIT_VER}},${GIT_VER},g
}

if [ -e src/version.h ]; then
	if (diff <(update_ver) src/version.h >/dev/null) then
		exit 0;
	fi
fi

echo "Updating 'src/version.h' ..." >/dev/stderr
update_ver > src/version.h
exit 0;
