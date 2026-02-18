import math

def is_valid(gridlines_i, totalline, totaldomain, error):
    """Validate the gridlines in a specific direction."""
    lnum = len(gridlines_i)

    if lnum == 0:
        error.append('[Error] No gridlines provided for validation.')
        return False

    if not math.isclose((gridlines_i[lnum - 1]['coord_f'] - gridlines_i[0]['coord_i']), totaldomain, rel_tol=1e-9):
        error.append(f"[Error] Computational domain length in {gridlines_i[0]['dir'].lower()}-direction does not match with the current gridline inputs.")
    if gridlines_i[0]['index_i'] != 1:
        error.append(f"[Error] The first index number is not equal to 1 in {gridlines_i[0]['dir'].lower()}-direction.")
    if gridlines_i[lnum - 1]['index_f'] != totalline:
        error.append(f"[Error] The final index number does not match the total gridline number in {gridlines_i[0]['dir'].lower()}-direction.")

    for n in range(1, lnum):
        if not math.isclose(gridlines_i[n]['coord_i'], gridlines_i[n - 1]['coord_f'], rel_tol=1e-9):
            error.append(f"[Error] The initial coordinate value of line {n + 1} in {gridlines_i[0]['dir'].lower()}-direction does not match the final coordinate value of the previous line.")
        if gridlines_i[n]['index_i'] != gridlines_i[n - 1]['index_f']:
            error.append(f"[Error] The initial index of line {n + 1} in {gridlines_i[0]['dir'].lower()}-direction does not match the final index of the previous line.")

    return not error