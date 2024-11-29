DEFAULT_SETTINGS = {
                    "Movies": True,
                    "WriteMoviesEvery": 1000,
                    }


class Generate:
    def __init__(self, **kwargs):
        config = DEFAULT_SETTINGS.copy()
        config.update(kwargs)
        for key, value in config.items():
            setattr(self, key, value)

    def WriteMovies(self):
        if self.Movies:
            return f"Movies\tyes\nWriteMoviesEvery\t{self.WriteMoviesEvery}\n\n"
        else:
            return "Movies\tno\n\n"

    def generate_config(self):
        return self.WriteMovies()