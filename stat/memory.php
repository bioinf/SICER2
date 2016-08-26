#!/usr/bin/php
# Usage:
# ./memory.php ibezdvornykh
<?php
$data = array(); 
$x = 0;
$file = 'test.' . time() . '.json';
$user = @$_SERVER['argv'][1];
if (!@$user) exit;
echo "# Datafile: $file\n";
while (true) {
	$pids = array();
	exec("ps -u {$user} | grep \"python\" | awk {'print $1'}", $pids);
	foreach ($pids as $pid) {
		if (!@$data[$pid]) $data[$pid] = array();
		$mem = array();
		exec("cat /proc/$pid/status | grep VmSize | awk {'print $2,$3'}", $mem);
		$sum = 0;
		foreach ($mem as $m) $sum += (int) $m;

		$data[$pid][] = array(time(), $sum);
		echo "\r$pid\tMemory: $sum [{$x}] ";
	}
	usleep(500000);
	$x++;
	if ($x%10 == 0) {
		$f = fopen($file, 'w+');
		fputs($f, json_encode($data));
		fclose($f);
	}
}


