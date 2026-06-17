"""Debug-output writers (null-object pattern)."""


class NullDebugWriter:
    """Null object that does nothing but implements the DebugWriter interface."""

    def populate(self, label, field_name, data):
        """No-op populate method."""
        pass

    def write(self, label):
        """No-op write method."""
        pass

    def populate_flux_components(self, label, flux, elem):
        """No-op populate flux components method."""
        pass


class DebugWriter:
    """Enhanced debug writer with populate functionality."""

    def __init__(self, debug, writer, mem):
        self.debug = debug
        self.writer = writer
        self.mem = mem
        self._pending_populations = {}

    def populate(self, label, field_name, data):
        """Populate field data for a given label."""
        if self.debug and self.writer is not None:
            if label not in self._pending_populations:
                self._pending_populations[label] = []
            self._pending_populations[label].append((field_name, data))

    def write(self, label):
        """Write all data including any pending populations."""
        if self.debug and self.writer is not None:
            # Apply any pending populations for this label
            if label in self._pending_populations:
                for field_name, data in self._pending_populations[label]:
                    self.writer.populate(label, field_name, data)
                del self._pending_populations[label]

            # Write the data
            self.writer.write_all(self.mem, label)

    def populate_flux_components(self, label, flux, elem):
        """Helper method to populate flux components."""
        if self.debug and self.writer is not None:
            self.populate(label, "rhoYu", flux[0].rhoY)
            self.populate(label, "rhoYv", flux[1].rhoY)
            if elem.ndim == 3:
                self.populate(label, "rhoYw", flux[2].rhoY)


def create_debug_writer(debug, writer, mem):
    """Factory function to create appropriate debug writer."""
    if debug and writer is not None:
        return DebugWriter(debug, writer, mem)
    else:
        return NullDebugWriter()
