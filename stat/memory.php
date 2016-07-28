#!/usr/bin/php
# Usage:
# ./memory.php ibezdvornykh

<?php
$data = []; 
$x = 0;
$file = 'memtest.' . time() . '.json';
while (true) {
	$user = $_SERVER['argv'][1];
	$pids = [];
	exec("ps -u {$user} | grep \"python\" | awk {'print $1'}", $pids);
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
		$f = fopen($file, 'w+');
		fputs($f, json_encode($data));
		fclose($f);
		echo "> $file; Saved!\n";
	}
}

