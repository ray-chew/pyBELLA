import h5py
import re

def do(h5path: str, ud, output_intervals: bool = False, time_increment: bool = False) -> None:
    if ud.diag_updt_targets:
        stepmax = ud.stepmax

        if output_intervals:
            # Determine steps to keep
            if stepmax <= 101:
                keep_steps = list(range(10, stepmax, 10))
            elif 101 < stepmax <= 1000:
                keep_steps = list(range(100, stepmax, 100))
            else:
                raise ValueError(f"stepmax {stepmax} > 1000 is too large as a diagnostic target.")
        else:
            keep_steps = [stepmax - 1]

        # If time_increment is True, add the step before each kept step
        if time_increment:
            additional_steps = [step - 1 for step in keep_steps if step > 0]
            keep_steps.extend(additional_steps)
            keep_steps = sorted(list(set(keep_steps)))

        output_path = h5path.replace(".h5", "_stripped.h5")

        with h5py.File(h5path, "r") as src, h5py.File(output_path, "w") as dst:
            # Copy all top-level groups, filtering time-tagged datasets
            for group in src.keys():
                src_group = src[group]
                dst_group = dst.create_group(group)

                for dset_name in src_group.keys():
                    # Identify timestep in dataset name, e.g. p2_nodes_010_after_full_step
                    step_match = re.search(r'_(\d+)_after_full_step$', dset_name)
                    if step_match:
                        step = int(step_match.group(1))
                        if step in keep_steps:
                            src_group.copy(dset_name, dst_group)
                    else:
                        # Copy non-timestep fields like "p2_nodes_ic"
                        src_group.copy(dset_name, dst_group)

        print(f"Stripped target file written to: {output_path}")