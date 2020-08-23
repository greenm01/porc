"""Just a simple pytest.

A sample test to get things started and to demonstrate proper pydoc strings.
Use google recommended style: https://google.github.io/styleguide/pyguide.html
"""


def return_it(obj: object):
    """Returns the argument passed in.

    Does nothing but simply return whatever is passed in.

    # Parameters:
    #    obj (object):Any object you want to pass.

    # Returns:
    #    Simply return the object that was passed in.
    """

    return obj


def test_foo():
    """A simple test that always passes.

    Just to make sure pytest is passing.
    """

    test_object = {}
    assert return_it(test_object) == test_object
