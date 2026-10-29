import math
from pathlib import Path
import selectors
import subprocess
import sys


driver = Path(sys.argv[1]).resolve()
build = Path(sys.argv[2]).resolve()
fixture = build / "ideal_rad.cache"
original = fixture.read_text().splitlines()
checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


def run(cache=fixture, query="", composition=(".7", ".02"), sha=None):
    command = [str(driver), str(cache), *composition]
    if sha is not None:
        command.append(sha)
    return subprocess.run(command, input=query, text=True, capture_output=True, timeout=10)


def reject(lines, label, status=11):
    filename = build / "malformed.cache"
    filename.write_text("\n".join(lines) + "\n")
    result = run(filename)
    require(result.returncode == 1, label + " exit")
    require(f"status={status}" in result.stderr, label + " status")
    require(not result.stdout, label + " no adopted table")


def replaced(line, token, value):
    lines = original.copy()
    tokens = lines[line].split()
    tokens[token] = value
    lines[line] = " ".join(tokens)
    return lines


for length in range(len(original)):
    reject(original[:length], f"truncated after {length} lines")
reject(["legacy"], "legacy cache")
reject(["FREEEOS_BIQUINTIC_V2", *original[1:]], "unknown version")
reject(original[:4] + original[5:], "old V1 missing provenance")
reject(original + ["extra"], "trailing data")
reject(original + [""], "trailing empty record")
reject(["x" * 4097, *original[1:]], "oversized record")
reject([" " * 400 + "x" * 129, *original[1:]], "oversized token")
for token in ("1", "0", "-1", "100001", "99999999999999999999", "2.0", "2*3"):
    reject(replaced(1, 0, token), "invalid dimension " + token)
reject(replaced(1, 0, "100000"), "bounded total nodes")
for line in (2, 5, 6):
    numeric_index = 3 if line == 6 else 0
    for value in ("NaN", "Inf", "1e9999", "1/", "2*1", "1,2", "1!comment"):
        reject(replaced(line, numeric_index, value), f"invalid number {line}:{value}")
for line in (2, 5):
    for token in range(len(original[line].split())):
        reject(replaced(line, token, ""), f"missing header token {line}:{token}")
        reject(replaced(line, token, original[line].split()[token] + " extra"), f"extra header token {line}:{token}")
for line, token, value in ((2, 0, "-0.1"), (2, 1, "-0.1"), (2, 0, ".99"), (2, 2, "-1"),
                            (2, 3, "2"), (2, 3, "1.0"), (5, 1, "0"), (5, 3, "-1"),
                            (5, 0, "1e300"), (5, 1, "1e-300"), (5, 3, "1e308"),
                            (4, 1, "0"), (4, 1, "4"), (4, 1, "2.0"), (4, 2, "MEAN"),
                            (4, 2, "unknown"), (4, 2, "mean,first")):
    reject(replaced(line, token, value), f"invalid header {line}:{token}:{value}")
for line in (3, 4):
    for sha in ("a" * 63, "a" * 65, "g" * 64, "a" * 63 + "/"):
        reject(replaced(line, 0, sha), "invalid SHA " + str(line))
for line, token, value in ((6, 0, "2"), (6, 1, "0"), (6, 2, "2"), (6, 2, "1.0"), (6, 2, "0")):
    reject(replaced(line, token, value), f"invalid node {token}:{value}")
for derivative in range(3, 12):
    reject(replaced(6, derivative, "NaN"), "nonfinite node derivative " + str(derivative))
for line_index in range(6, len(original)):
    reject(original[:line_index] + original[line_index + 1:], "missing node " + str(line_index))
swapped = original.copy()
swapped[6], swapped[7] = swapped[7], swapped[6]
reject(swapped, "out of order nodes")
for radius in ("1", "2", "3"):
    for route in ("mean", "first", "second"):
        lines = replaced(4, 1, radius)
        tokens = lines[4].split()
        tokens[2] = route
        lines[4] = " ".join(tokens)
        filename = build / "valid_processing.cache"
        filename.write_text("\n".join(lines) + "\n")
        require(run(filename).returncode == 0, "accepted provenance " + radius + route)
uppercase = original.copy()
uppercase[3] = uppercase[3].upper()
uppercase[4] = uppercase[4].upper().replace("MEAN", "mean")
filename = build / "uppercase_sha.cache"
filename.write_text("\n".join(uppercase) + "\n")
require(run(filename, sha="A" * 64).returncode == 0, "case-insensitive SHA")
require(run(composition=(".6", ".02")).returncode == 1, "composition mismatch")
require(run(sha="b" * 64).returncode == 1, "source mismatch")
require("status=14" in run(sha="bad").stderr, "invalid expected hash")
require("status=10" in run(build / "does_not_exist.cache").stderr, "missing cache")

offset = 12 * math.log(10)
temperature_coordinate = 5.05
query_x = 0.2 + 2 * temperature_coordinate - offset
query = f"{query_x:.17g} {temperature_coordinate:.17g}\n"
result = run(query=query * 256 + "-100 5.05\n")
require(result.returncode == 0, "batch process")
records = result.stdout.splitlines()
names = records[0].split()
require(len(names) == 58 and len(set(names)) == 58, "stable 58-column unique header")
require(names[:5] == ["x", "t", "status_real", "status_ad", "status_state"], "header statuses")
require(len(records) == 258, "one row per query")
for record in records[1:-1]:
    values = list(map(float, record.split()))
    require(len(values) == 58 and values[2:5] == [0, 0, 0], "successful batch status")
    require(all(math.isfinite(value) for value in values), "finite batch output")
rejected = list(map(float, records[-1].split()))
require(rejected[2:5] == [15, 15, 15], "out of bounds status row")
require(all(math.isnan(value) for value in rejected[5:]), "failed row NaN poisoning")
for invalid_query in ("", "1\n", "1 2 3\n", "1*2 5\n", "1/ 5\n", "1,2 5\n", "NaN 5\n", " " * 4097 + "\n"):
    result = run(query=invalid_query)
    require(result.returncode == (0 if not invalid_query else 2), "strict query record " + repr(invalid_query[:20]))
require(run(query=query.rstrip()).returncode == 0, "final record without newline")
result = run(query="1e9999 5\n")
require(result.returncode == 0 and result.stdout.splitlines()[1].split()[2:5] == ["14", "14", "14"], "nonfinite query status")

for native_u in (0.0, 0.4, 0.8):
    for native_t in (5.0, 5.1):
        physical_x = native_u + 2 * native_t - offset
        result = run(query=f"{physical_x:.17g} {native_t:.17g}\n")
        require(result.stdout.splitlines()[1].split()[2:5] == ["0", "0", "0"], "physical endpoint round trip")
for native_u in (-1e-10, 0.8 + 1e-10):
    physical_x = native_u + 2 * temperature_coordinate - offset
    result = run(query=f"{physical_x:.17g} {temperature_coordinate:.17g}\n")
    require(result.stdout.splitlines()[1].split()[2:5] == ["15", "15", "15"], "outside endpoint envelope rejected")
for native_t in (math.nextafter(5.0, -math.inf), math.nextafter(5.1, math.inf)):
    physical_x = 0.2 + 2 * native_t - offset
    result = run(query=f"{physical_x:.17g} {native_t:.17g}\n")
    require(result.stdout.splitlines()[1].split()[2:5] == ["15", "15", "15"], "temperature boundary strictly enforced")

with subprocess.Popen([str(driver), str(fixture), ".7", ".02"], stdin=subprocess.PIPE,
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1) as process:
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    require(bool(selector.select(5)), "header flushed before input")
    require(process.stdout.readline().split() == names, "persistent header")
    for batch in range(3):
        process.stdin.write(query)
        process.stdin.flush()
        require(bool(selector.select(5)), "query output flushed before EOF")
        require(len(process.stdout.readline().split()) == 58, "persistent numerical row")
    process.stdin.close()
    require(process.wait(timeout=5) == 0, "persistent EOF")
    selector.close()

print(f"PASS {checks} standalone driver/cache protocol checks")
