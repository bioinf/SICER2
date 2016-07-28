#!/usr/bin/php
<?php
# ./mem.php
$data = [];
$x = 0;
while (true) {
	$pids = [];
	exec("ps -u bezdvornyhi | grep \"python\" | awk {'print $1'}", $pids);
	foreach ($pids as $pid) {
		if (!$data[$pid]) $data[$pid] = [];
		$mem = []; $sum = 0;
		exec("cat /proc/$pid/status | grep VmSize | awk {'print $2,$3'}", $mem);
		foreach ($mem as $m) $sum += (int) $m;
		$data[$pid][] = [time(), $m];
		echo "$pid\tMemory: $m\n";
	}
	usleep(100000);
	$x++;
	if ($x%10 == 0) {
		$f = fopen('memory.txt', 'w+');
		fputs($f, json_encode($data));
		fclose($f);
		echo "Saved!\n";
	}
}

