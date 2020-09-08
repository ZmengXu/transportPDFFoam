# GRI 2.11

For the Methane flames in the tutorials the GRI 2.11 mechanism was used. This
is not the current versions, which is 3.0, however for the purpose of these
flames, the impact should be minimal. The reason this older version was used
is historical, as in the work of B. Zoller also NOx-formation was an issue,
and in that field, version 2.11 performs significantly better than the newer
3.0 version.

The mechanism was obtained from the [GRI website](http://diesel.me.berkeley.edu/~gri_mech/new21/version21/text21.html)
It has then been converted to the FlameMaster format using the tools available
in the FlameMaster distribution, available from H. Pitsch [[FlameMaster]](../../../doc/references.md#FlameMaster).
The process is as follows

```sh
$ CreateBinFile -i thermo211.dat -m transport.dat -o gri.211.bin
$ ScanMan -i grimech211.dat -t gri.211.bin -o gri.211.pre -S
```
