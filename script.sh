
echo bench-l1-g8 >> log.txt
time mpirun -n 22 bin/athena -i bench-l1-g8   >> log.txt
echo bench-l1-g16 >> log.txt
time mpirun -n 22 bin/athena -i bench-l1-g16  >> log.txt
echo bench-l1-g32 >> log.txt
time mpirun -n 22 bin/athena -i bench-l1-g32  >> log.txt
echo bench-l1-g64 >> log.txt
time mpirun -n 22 bin/athena -i bench-l1-g64  >> log.txt

echo bench-l3-g8 >> log.txt
time mpirun -n 22 bin/athena -i bench-l3-g8   >> log.txt
echo bench-l3-g16 >> log.txt
time mpirun -n 22 bin/athena -i bench-l3-g16  >> log.txt
echo bench-l3-g32 >> log.txt
time mpirun -n 22 bin/athena -i bench-l3-g32  >> log.txt
echo bench-l3-g64 >> log.txt
time mpirun -n 22 bin/athena -i bench-l3-g64  >> log.txt

echo bench-l6-g8 >> log.txt
time mpirun -n 22 bin/athena -i bench-l6-g8   >> log.txt
echo bench-l6-g16 >> log.txt
time mpirun -n 22 bin/athena -i bench-l6-g16  >> log.txt
echo bench-l6-g32 >> log.txt
time mpirun -n 22 bin/athena -i bench-l6-g32  >> log.txt
echo bench-l6-g64 >> log.txt
time mpirun -n 22 bin/athena -i bench-l6-g64  >> log.txt
