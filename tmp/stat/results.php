#!/usr/bin/php
<?php
$data = [];
$json = json_decode(file_get_contents('test.1472265719.json'));

exec("ls -l /media/DISK1/SICER_proj/hs_marks/bams/*.bam | awk '{print $5,$9}'", $fx);
foreach ($fx as $l) {
	list($size, $filename) = explode(' ', $l);
	$filename = explode('/', $filename);
	$filename = end($filename);
	$filename = str_replace('.bam', '', $filename);
	if ($json->{$filename}) {
		$list = [];
		foreach ($json->{$filename} as $tx) $list[] = $tx[1];
		$data[] = [max($list)/2024, $size/2024/2024];
	}
}
echo json_encode($data);

exit;
$data = [];
exec("ls -l /media/DISK1/SICER_proj/hs_marks/bams/*.bam | awk '{print $5,$9}'", $fx);
foreach ($fx as $l) {
	list($size, $filename) = explode(' ', $l);
	$filename = explode('/', $filename);
	$filename = end($filename);
	$filename = '/media/DISK1/SICER_proj/hs_marks/new_SICER_x1/logs/' . str_replace('.bam', '.200.log', $filename);
	if (!file_exists($filename)) continue;
	$time = [];
	exec("cat $filename | grep \"seconds\"", $time);
	if (@$time[0]) $data[] = [(float) explode(': ', $time[0])[1], $size/2024/2024];
}
echo json_encode($data);

